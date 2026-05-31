//! Shannon's noisy channel coding theorem bounds.


/// Binary entropy function H(p) = -p*log2(p) - (1-p)*log2(1-p).
pub fn binary_entropy(p: f64) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        return 0.0;
    }
    -p * p.log2() - (1.0 - p) * (1.0 - p).log2()
}

/// Channel capacity of a Binary Symmetric Channel (BSC) with crossover probability p.
/// C = 1 - H(p) bits per channel use.
pub fn bsc_capacity(p: f64) -> f64 {
    1.0 - binary_entropy(p)
}

/// Shannon limit: maximum code rate R for reliable communication over BSC(p).
/// R < C = 1 - H(p).
pub fn shannon_limit_bsc(p: f64) -> f64 {
    bsc_capacity(p)
}

/// Minimum Eb/N0 (in dB) required for reliable communication at rate R over BSC(p).
/// For BSC: Eb/N0 = (1/(2*R)) * ((1-2p)^(-2))... 
/// This is a simplified approximation.
pub fn min_eb_n0_db(rate: f64, p: f64) -> f64 {
    if p <= 0.0 || p >= 0.5 || rate <= 0.0 {
        return f64::INFINITY;
    }
    // For BSC: capacity = 1 - H(p)
    // Shannon limit for BPSK: Eb/N0 = (2^R - 1) approximation
    // More precisely, for BSC: Eb/N0 ≈ (1/(2*R)) * ln(1/(1-2p)) approximation
    let eb_n0 = 1.0 / (2.0 * rate) * (1.0 / (1.0 - 2.0 * p)).ln();
    10.0 * eb_n0.log10()
}

/// Gilbert-Varshamov bound: maximum code length n for given d and k.
/// A linear [n, k] code with minimum distance at least d exists if:
/// sum_{i=0}^{d-2} C(n-1, i) < 2^(n-k)
pub fn gilbert_varshamov_bound(n: usize, d: usize, k: usize) -> bool {
    let mut sum: f64 = 0.0;
    for i in 0..=(d - 2) {
        sum += binomial(n - 1, i) as f64;
    }
    let rhs = 2.0_f64.powi((n - k) as i32);
    sum < rhs
}

/// Hamming bound (sphere-packing bound): A q-ary (n, M, d) code satisfies:
/// M * sum_{i=0}^{t} C(n, i) * (q-1)^i <= q^n, where t = floor((d-1)/2).
pub fn hamming_bound(n: usize, d: usize, q: usize, num_codewords: usize) -> bool {
    let t = (d - 1) / 2;
    let qf = q as f64;
    let mut volume: f64 = 0.0;
    for i in 0..=t {
        volume += binomial(n, i) as f64 * (qf - 1.0).powi(i as i32);
    }
    (num_codewords as f64) * volume <= qf.powi(n as i32)
}

/// Singleton bound: d_min <= n - k + 1.
pub fn singleton_bound(n: usize, k: usize) -> usize {
    n - k + 1
}

/// Check if a code is MDS (Maximum Distance Separable): d_min = n - k + 1.
pub fn is_mds(n: usize, k: usize, d_min: usize) -> bool {
    d_min == singleton_bound(n, k)
}

/// Compute binomial coefficient C(n, k).
pub fn binomial(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result: u64 = 1;
    for i in 0..k {
        result = result * (n - i) as u64 / (i + 1) as u64;
    }
    result
}

/// Sphere-packing bound for binary codes: maximum number of codewords.
/// M <= 2^n / sum_{i=0}^{t} C(n, i)
pub fn max_codewords_sphere_packing(n: usize, t: usize) -> f64 {
    let mut volume: f64 = 0.0;
    for i in 0..=t {
        volume += binomial(n, i) as f64;
    }
    (1u64 << n) as f64 / volume
}

/// Plotkin bound: for binary codes with d > n/2:
/// M <= 2d / (2d - n).
pub fn plotkin_bound(n: usize, d: usize) -> f64 {
    if 2 * d <= n {
        return f64::INFINITY; // Bound doesn't apply
    }
    (2 * d) as f64 / (2 * d - n) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_entropy() {
        assert_eq!(binary_entropy(0.0), 0.0);
        assert_eq!(binary_entropy(1.0), 0.0);
        assert!((binary_entropy(0.5) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_binary_entropy_mid() {
        let h = binary_entropy(0.25);
        assert!(h > 0.0 && h < 1.0);
        assert!((h - 0.81127812445).abs() < 1e-8);
    }

    #[test]
    fn test_bsc_capacity() {
        assert_eq!(bsc_capacity(0.0), 1.0); // Perfect channel
        assert!((bsc_capacity(0.5) - 0.0).abs() < 1e-10); // Useless channel
    }

    #[test]
    fn test_shannon_limit() {
        let limit = shannon_limit_bsc(0.1);
        assert!(limit > 0.0 && limit < 1.0);
    }

    #[test]
    fn test_singleton_bound() {
        assert_eq!(singleton_bound(7, 4), 4);
        assert_eq!(singleton_bound(15, 11), 5);
    }

    #[test]
    fn test_is_mds() {
        // Reed-Solomon codes are MDS
        assert!(is_mds(7, 3, 5));
        // Hamming [7,4,3] is not MDS (Singleton gives 4)
        assert!(!is_mds(7, 4, 3));
    }

    #[test]
    fn test_hamming_bound() {
        // [7,4,3] Hamming code: M=16, t=1
        // 16 * (1 + 7) = 128 = 2^7 -> tight!
        assert!(hamming_bound(7, 3, 2, 16));
    }

    #[test]
    fn test_gilbert_varshamov() {
        // [7,4] code with d=3 should satisfy GV bound
        assert!(gilbert_varshamov_bound(7, 3, 4));
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 1), 5);
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(10, 5), 252);
    }

    #[test]
    fn test_max_codewords_sphere_packing() {
        // For n=7, t=1: max = 128 / 8 = 16
        let max = max_codewords_sphere_packing(7, 1);
        assert!((max - 16.0).abs() < 1e-6);
    }

    #[test]
    fn test_plotkin_bound() {
        // n=7, d=4: 2d/(2d-n) = 8/1 = 8
        let pb = plotkin_bound(7, 4);
        assert!((pb - 8.0).abs() < 1e-6);
    }

    #[test]
    fn test_bsc_capacity_monotone() {
        let c1 = bsc_capacity(0.05);
        let c2 = bsc_capacity(0.1);
        let c3 = bsc_capacity(0.2);
        assert!(c1 > c2);
        assert!(c2 > c3);
    }

    #[test]
    fn test_hamming_bound_perfect_code() {
        // Perfect code: [23, 12, 7] Golay code, t=3
        // M = 2^12 = 4096
        // Volume = 1 + 23 + 253 + 1771 = 2048
        // 4096 * 2048 = 2^12 * 2^11 = 2^23 ✓
        assert!(hamming_bound(23, 7, 2, 4096));
    }
}
