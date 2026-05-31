//! Reed-Solomon basics: encoding and syndrome computation over GF(2^m).

use crate::GFElement;

/// Reed-Solomon code parameters.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ReedSolomonCode {
    /// Field parameter m (GF(2^m)).
    pub m: u8,
    /// Codeword length n.
    pub n: usize,
    /// Data length k.
    pub k: usize,
    /// Number of parity symbols t = (n - k) / 2.
    pub t: usize,
    /// First consecutive root (usually 1).
    pub first_root: usize,
}

impl ReedSolomonCode {
    /// Create a Reed-Solomon code over GF(2^m) with n codeword symbols and k data symbols.
    pub fn new(m: u8, n: usize, k: usize) -> Self {
        let t = (n - k) / 2;
        ReedSolomonCode { m, n, k, t, first_root: 1 }
    }

    /// Create a Reed-Solomon code with custom first root.
    pub fn with_first_root(m: u8, n: usize, k: usize, first_root: usize) -> Self {
        let t = (n - k) / 2;
        ReedSolomonCode { m, n, k, t, first_root }
    }

    /// Get primitive element alpha (value 2 in GF(2^m)).
    pub fn alpha(&self) -> GFElement {
        GFElement::new(2, self.m)
    }

    /// Encode data symbols into codeword using systematic encoding.
    /// c(x) = d(x) * x^(n-k) + r(x) where r(x) = d(x) * x^(n-k) mod g(x).
    pub fn encode(&self, data: &[GFElement]) -> Vec<GFElement> {
        assert_eq!(data.len(), self.k);
        let parity_count = self.n - self.k;
        let alpha = self.alpha();
        let zero = GFElement::new(0, self.m);

        // Build generator polynomial g(x) = prod_{i=first_root}^{first_root+2t-1} (x - alpha^i)
        let mut gen = vec![GFElement::new(1, self.m)]; // Start with 1
        for j in 0..2 * self.t {
            let root = alpha.pow((self.first_root + j) as u32);
            // Multiply gen by (x - root) = (x + root) in GF(2^m)
            let mut new_gen = vec![zero; gen.len() + 1];
            for i in 0..gen.len() {
                new_gen[i] = new_gen[i].add(&gen[i].mul(&root)); // shift-multiply by root
                new_gen[i + 1] = new_gen[i + 1].add(&gen[i]);  // shift (x)
            }
            gen = new_gen;
        }

        // Compute d(x) * x^(n-k)
        let mut shifted = vec![zero; parity_count];
        shifted.extend_from_slice(data);

        // Polynomial division: shifted / gen -> remainder
        let gen_deg = gen.len() - 1;
        let mut remainder = shifted.clone();
        for i in (gen_deg..remainder.len()).rev() {
            if remainder[i].val != 0 {
                let coeff = remainder[i].div(&gen[gen_deg]).unwrap();
                for j in 0..=gen_deg {
                    remainder[i - gen_deg + j] = remainder[i - gen_deg + j].sub(&gen[j].mul(&coeff));
                }
            }
        }

        // Parity is the remainder (lowest parity_count coefficients)
        let mut codeword = data.to_vec();
        for i in 0..parity_count {
            codeword.push(remainder[i]);
        }
        codeword
    }

    /// Compute syndromes S_j = received(alpha^(first_root + j)) for j = 0..2t.
    pub fn syndromes(&self, received: &[GFElement]) -> Vec<GFElement> {
        let alpha = self.alpha();
        let num_syndromes = 2 * self.t;
        let mut syndromes = Vec::with_capacity(num_syndromes);

        for j in 0..num_syndromes {
            let alpha_j = alpha.pow((self.first_root + j) as u32);
            // Evaluate received polynomial at alpha_j using Horner's method
            let mut eval = GFElement::new(0, self.m);
            for i in (0..received.len()).rev() {
                eval = alpha_j.mul(&eval).add(&received[i]);
            }
            syndromes.push(eval);
        }

        syndromes
    }

    /// Check if syndromes indicate no errors.
    pub fn is_valid(&self, received: &[GFElement]) -> bool {
        self.syndromes(received).iter().all(|s| s.val == 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gf(val: u16) -> GFElement {
        GFElement::new(val, 4)
    }

    #[test]
    fn test_rs_syndromes_no_error() {
        // RS(15, 11) over GF(2^4), t = 2
        let rs = ReedSolomonCode::new(4, 15, 11);
        // All-zeros codeword should have zero syndromes
        let codeword = vec![gf(0); 15];
        let syndromes = rs.syndromes(&codeword);
        assert!(syndromes.iter().all(|s| s.val == 0));
    }

    #[test]
    fn test_rs_encode_syndromes() {
        let rs = ReedSolomonCode::new(4, 15, 11);
        let data: Vec<GFElement> = (0..11).map(|i| gf(i)).collect();
        let codeword = rs.encode(&data);
        assert_eq!(codeword.len(), 15);
        let syndromes = rs.syndromes(&codeword);
        assert!(syndromes.iter().all(|s| s.val == 0), "Syndromes should be zero for valid codeword");
    }

    #[test]
    fn test_rs_detects_error() {
        let rs = ReedSolomonCode::new(4, 15, 11);
        let data: Vec<GFElement> = (0..11).map(|i| gf(i)).collect();
        let mut codeword = rs.encode(&data);

        // Introduce an error
        codeword[3] = codeword[3].add(&gf(1));

        let syndromes = rs.syndromes(&codeword);
        assert!(syndromes.iter().any(|s| s.val != 0), "Syndromes should be nonzero for corrupted codeword");
    }

    #[test]
    fn test_rs_parameters() {
        let rs = ReedSolomonCode::new(4, 15, 11);
        assert_eq!(rs.t, 2);
        assert_eq!(rs.n - rs.k, 4); // 2t parity symbols
    }

    #[test]
    fn test_rs_gf3_basic() {
        let rs = ReedSolomonCode::new(3, 7, 3);
        assert_eq!(rs.t, 2);
        let data = vec![GFElement::new(0, 3), GFElement::new(1, 3), GFElement::new(2, 3)];
        let codeword = rs.encode(&data);
        assert_eq!(codeword.len(), 7);
        let syndromes = rs.syndromes(&codeword);
        assert!(syndromes.iter().all(|s| s.val == 0));
    }

    #[test]
    fn test_rs_gf3_error_detection() {
        let rs = ReedSolomonCode::new(3, 7, 3);
        let data = vec![GFElement::new(0, 3), GFElement::new(1, 3), GFElement::new(2, 3)];
        let mut codeword = rs.encode(&data);
        codeword[0] = codeword[0].add(&GFElement::new(1, 3));
        assert!(!rs.is_valid(&codeword));
    }

    #[test]
    fn test_rs_all_zeros() {
        let rs = ReedSolomonCode::new(4, 15, 11);
        let data = vec![gf(0); 11];
        let codeword = rs.encode(&data);
        assert!(codeword.iter().all(|c| c.val == 0));
    }
}
