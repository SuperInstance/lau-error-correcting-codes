//! Linear codes: generator matrices, parity check, minimum distance.

use crate::{Bit, BitVec, gf2};

/// A linear code defined by a generator matrix.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct LinearCode {
    /// Generator matrix G (k x n): rows are basis vectors.
    pub generator: Vec<BitVec>,
    /// Codeword length n.
    pub n: usize,
    /// Data length k.
    pub k: usize,
}

impl LinearCode {
    /// Create a linear code from a generator matrix.
    pub fn new(generator: Vec<BitVec>) -> Self {
        let k = generator.len();
        let n = generator.first().map(|r| r.len()).unwrap_or(0);
        LinearCode { generator, n, k }
    }

    /// Create a linear code from a parity check matrix H (r x n).
    /// Finds a generator G such that G * H^T = 0.
    pub fn from_parity_check(h: Vec<BitVec>) -> Self {
        let r = h.len();
        let n = h.first().map(|row| row.len()).unwrap_or(0);
        let _k = n - r;

        // Build systematic generator: G = [I_k | P] where H = [P^T | I_r]
        // We need to find the systematic form of H, then extract G.
        // For simplicity, we'll build G from H using Gaussian elimination.

        // Augment H and perform row reduction to systematic form
        let mut h_mat = h.clone();

        // Gaussian elimination on H
        let mut pivot_cols: Vec<usize> = Vec::new();
        for row in 0..r {
            // Find pivot column
            let mut found = false;
            for col in 0..n {
                if h_mat[row][col] == 1 && !pivot_cols.contains(&col) {
                    pivot_cols.push(col);
                    // Eliminate this column from other rows
                    for other in 0..r {
                        if other != row && h_mat[other][col] == 1 {
                            for c in 0..n {
                                h_mat[other][c] ^= h_mat[row][c];
                            }
                        }
                    }
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }

        // Data columns are non-pivot columns
        let data_cols: Vec<usize> = (0..n).filter(|c| !pivot_cols.contains(c)).collect();

        // Build generator matrix
        let mut g = Vec::new();
        for &dc in &data_cols {
            let mut row = vec![0u8; n];
            row[dc] = 1;
            for (i, &pc) in pivot_cols.iter().enumerate() {
                row[pc] = h_mat[i][dc];
            }
            g.push(row);
        }

        LinearCode::new(g)
    }

    /// Encode data bits using the generator matrix.
    pub fn encode(&self, data: &[Bit]) -> BitVec {
        assert_eq!(data.len(), self.k, "Data must be {} bits", self.k);
        let mut codeword = vec![0u8; self.n];
        for i in 0..self.k {
            if data[i] == 1 {
                for j in 0..self.n {
                    codeword[j] ^= self.generator[i][j];
                }
            }
        }
        codeword
    }

    /// Compute the parity check matrix from the generator.
    /// Assumes G is in systematic form [I_k | P].
    pub fn parity_check_matrix(&self) -> Vec<BitVec> {
        let r = self.n - self.k;
        let mut h = Vec::with_capacity(r);

        // If systematic, G = [I_k | P], then H = [P^T | I_r]
        for j in 0..r {
            let mut row = vec![0u8; self.n];
            // P^T part
            for i in 0..self.k {
                row[i] = self.generator[i][self.k + j];
            }
            // I_r part
            row[self.k + j] = 1;
            h.push(row);
        }
        h
    }

    /// Compute syndrome: H * received^T.
    pub fn syndrome(&self, received: &[Bit]) -> BitVec {
        let h = self.parity_check_matrix();
        h.iter().map(|row| gf2::dot(row, received)).collect()
    }

    /// Check if received word is a valid codeword.
    pub fn is_valid(&self, received: &[Bit]) -> bool {
        self.syndrome(received).iter().all(|&b| b == 0)
    }

    /// Compute the minimum distance (Hamming weight of the minimum-weight non-zero codeword).
    /// Brute force: enumerate all 2^k codewords.
    pub fn minimum_distance(&self) -> usize {
        if self.k == 0 {
            return 0;
        }

        let total = 1usize << self.k;
        let mut min_weight = self.n;

        for i in 1..total {
            let data: BitVec = (0..self.k).map(|j| ((i >> j) & 1) as u8).collect();
            let codeword = self.encode(&data);
            let weight = codeword.iter().filter(|&&b| b == 1).count();
            if weight < min_weight {
                min_weight = weight;
            }
        }
        min_weight
    }

    /// Compute the weight (number of non-zero elements) of a vector.
    pub fn weight(v: &[Bit]) -> usize {
        v.iter().filter(|&&b| b == 1).count()
    }

    /// Compute Hamming distance between two vectors.
    pub fn hamming_distance(a: &[Bit], b: &[Bit]) -> usize {
        a.iter().zip(b.iter()).filter(|(&x, &y)| x != y).count()
    }

    /// Get all codewords (only practical for small k).
    pub fn all_codewords(&self) -> Vec<BitVec> {
        let total = 1usize << self.k;
        (0..total)
            .map(|i| {
                let data: BitVec = (0..self.k).map(|j| ((i >> j) & 1) as u8).collect();
                self.encode(&data)
            })
            .collect()
    }

    /// The error-detection capability: d_min - 1.
    pub fn error_detection_capability(&self) -> usize {
        self.minimum_distance().saturating_sub(1)
    }

    /// The error-correction capability: floor((d_min - 1) / 2).
    pub fn error_correction_capability(&self) -> usize {
        (self.minimum_distance().saturating_sub(1)) / 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_code_encode() {
        // Simple [5, 2] code
        let g = vec![
            vec![1, 0, 1, 1, 0],
            vec![0, 1, 0, 1, 1],
        ];
        let code = LinearCode::new(g);
        assert_eq!(code.n, 5);
        assert_eq!(code.k, 2);

        let codeword = code.encode(&[0, 0]);
        assert_eq!(codeword, vec![0, 0, 0, 0, 0]);

        let codeword = code.encode(&[1, 0]);
        assert_eq!(codeword, vec![1, 0, 1, 1, 0]);
    }

    #[test]
    fn test_linear_code_minimum_distance() {
        // Repetition code [3, 1]: G = [[1,1,1]], d_min = 3
        let g = vec![vec![1, 1, 1]];
        let code = LinearCode::new(g);
        assert_eq!(code.minimum_distance(), 3);
    }

    #[test]
    fn test_linear_code_min_distance_simple() {
        let g = vec![
            vec![1, 0, 1],
            vec![0, 1, 1],
        ];
        let code = LinearCode::new(g);
        // Codewords: 000, 101, 011, 110 -> min weight = 2
        assert_eq!(code.minimum_distance(), 2);
    }

    #[test]
    fn test_parity_check_systematic() {
        let g = vec![
            vec![1, 0, 0, 1, 1],
            vec![0, 1, 0, 0, 1],
            vec![0, 0, 1, 1, 0],
        ];
        let code = LinearCode::new(g);
        let h = code.parity_check_matrix();
        assert_eq!(h.len(), 2);
        // H * G^T should be 0
        for grow in &code.generator {
            for hrow in &h {
                assert_eq!(gf2::dot(grow, hrow), 0);
            }
        }
    }

    #[test]
    fn test_syndrome_decoding() {
        let g = vec![
            vec![1, 0, 0, 1],
            vec![0, 1, 0, 1],
            vec![0, 0, 1, 1],
        ];
        let code = LinearCode::new(g);
        let codeword = code.encode(&[1, 0, 1]);
        assert!(code.is_valid(&codeword));

        let mut corrupted = codeword.clone();
        corrupted[2] ^= 1;
        assert!(!code.is_valid(&corrupted));
        let syndrome = code.syndrome(&corrupted);
        assert!(syndrome.iter().any(|&b| b == 1));
    }

    #[test]
    fn test_from_parity_check() {
        let h = vec![
            vec![1, 1, 0],
            vec![0, 1, 1],
        ];
        let code = LinearCode::from_parity_check(h);
        assert_eq!(code.n, 3);
        assert_eq!(code.k, 1);
        // The generator should produce valid codewords
        for i in 0..2 {
            let data = vec![i as u8];
            let cw = code.encode(&data);
            assert!(code.is_valid(&cw), "Codeword {:?} should be valid", cw);
        }
    }

    #[test]
    fn test_error_capabilities() {
        let g = vec![vec![1, 1, 1]];
        let code = LinearCode::new(g);
        assert_eq!(code.error_detection_capability(), 2);
        assert_eq!(code.error_correction_capability(), 1);
    }

    #[test]
    fn test_all_codewords() {
        let g = vec![
            vec![1, 0, 1],
            vec![0, 1, 1],
        ];
        let code = LinearCode::new(g);
        let cws = code.all_codewords();
        assert_eq!(cws.len(), 4);
    }

    #[test]
    fn test_hamming_distance() {
        assert_eq!(LinearCode::hamming_distance(&[0, 0, 0], &[0, 0, 0]), 0);
        assert_eq!(LinearCode::hamming_distance(&[1, 0, 1], &[0, 1, 1]), 2);
    }
}
