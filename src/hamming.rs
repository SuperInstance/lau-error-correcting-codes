//! Hamming codes: encoding, syndrome decoding, error correction.

use crate::{Bit, BitVec, gf2};

/// A Hamming code for a given number of parity bits r.
/// Codeword length n = 2^r - 1, data bits k = 2^r - r - 1.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct HammingCode {
    /// Number of parity bits.
    pub r: usize,
    /// Codeword length.
    pub n: usize,
    /// Data bits.
    pub k: usize,
}

impl HammingCode {
    /// Create a Hamming code with r parity bits.
    pub fn new(r: usize) -> Self {
        assert!(r >= 2, "Hamming codes require at least 2 parity bits");
        let n = (1 << r) - 1;
        let k = n - r;
        HammingCode { r, n, k }
    }

    /// Build the parity check matrix H (r x n).
    /// Column j contains the binary representation of (j+1).
    pub fn parity_check_matrix(&self) -> Vec<BitVec> {
        let mut h = Vec::new();
        for i in 0..self.r {
            let mut row = BitVec::with_capacity(self.n);
            for j in 1..=self.n {
                row.push(((j >> i) & 1) as Bit);
            }
            h.push(row);
        }
        h
    }

    /// Build the generator matrix G (k x n) in systematic form.
    pub fn generator_matrix(&self) -> Vec<BitVec> {
        let h = self.parity_check_matrix();
        // In systematic form: H = [P^T | I_r], G = [I_k | P]
        // Find which positions are parity (powers of 2) and data
        let parity_positions: Vec<usize> = (0..self.r).map(|i| (1 << i) - 1).collect();
        let data_positions: Vec<usize> = (0..self.n)
            .filter(|j| !parity_positions.contains(j))
            .collect();

        let mut g = Vec::new();
        for row_idx in 0..self.k {
            let mut row = vec![0u8; self.n];
            // Identity in data positions
            row[data_positions[row_idx]] = 1;
            // P part: for each parity position, compute the dot product
            for (pi, &pp) in parity_positions.iter().enumerate() {
                row[pp] = h[pi][data_positions[row_idx]];
            }
            g.push(row);
        }
        g
    }

    /// Encode k data bits into an n-bit codeword.
    pub fn encode(&self, data: &[Bit]) -> BitVec {
        assert_eq!(data.len(), self.k, "Data must be {} bits", self.k);
        let g = self.generator_matrix();
        let mut codeword = vec![0u8; self.n];
        for j in 0..self.n {
            for i in 0..self.k {
                codeword[j] ^= data[i] & g[i][j];
            }
        }
        codeword
    }

    /// Compute syndrome of received word.
    pub fn syndrome(&self, received: &[Bit]) -> BitVec {
        assert_eq!(received.len(), self.n, "Received must be {} bits", self.n);
        let h = self.parity_check_matrix();
        h.iter().map(|row| gf2::dot(row, received)).collect()
    }

    /// Decode using syndrome decoding. Returns corrected codeword and error position (None if no error).
    pub fn decode(&self, received: &[Bit]) -> (BitVec, Option<usize>) {
        let syndrome = self.syndrome(received);
        let syndrome_val: usize = syndrome.iter().enumerate().fold(0usize, |acc, (i, &b)| {
            acc | ((b as usize) << i)
        });

        let mut corrected = received.to_vec();
        let error_pos = if syndrome_val == 0 {
            None
        } else if syndrome_val <= self.n {
            corrected[syndrome_val - 1] ^= 1;
            Some(syndrome_val - 1)
        } else {
            // Uncorrectable (more than 1 error for standard Hamming)
            None
        };

        (corrected, error_pos)
    }

    /// Extract data bits from a corrected codeword.
    pub fn extract_data(&self, codeword: &[Bit]) -> BitVec {
        let parity_positions: Vec<usize> = (0..self.r).map(|i| (1 << i) - 1).collect();
        codeword.iter().enumerate()
            .filter(|(i, _)| !parity_positions.contains(i))
            .map(|(_, &b)| b)
            .collect()
    }

    /// Full decode: correct errors and return data bits.
    pub fn decode_to_data(&self, received: &[Bit]) -> BitVec {
        let (corrected, _) = self.decode(received);
        self.extract_data(&corrected)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_7_4_encode() {
        let code = HammingCode::new(3);
        assert_eq!(code.n, 7);
        assert_eq!(code.k, 4);
        let data = vec![1, 0, 1, 1];
        let codeword = code.encode(&data);
        // Verify it's a valid codeword
        let syndrome = code.syndrome(&codeword);
        assert!(syndrome.iter().all(|&b| b == 0));
    }

    #[test]
    fn test_hamming_7_4_no_error() {
        let code = HammingCode::new(3);
        let data = vec![1, 0, 1, 1];
        let codeword = code.encode(&data);
        let (corrected, pos) = code.decode(&codeword);
        assert_eq!(corrected, codeword);
        assert_eq!(pos, None);
    }

    #[test]
    fn test_hamming_7_4_single_error_correction() {
        let code = HammingCode::new(3);
        let data = vec![1, 0, 1, 1];
        let codeword = code.encode(&data);

        // Flip each bit and verify correction
        for i in 0..7 {
            let mut corrupted = codeword.clone();
            corrupted[i] ^= 1;
            let (corrected, pos) = code.decode(&corrupted);
            assert_eq!(corrected, codeword, "Failed to correct error at position {}", i);
            assert_eq!(pos, Some(i));
        }
    }

    #[test]
    fn test_hamming_7_4_decode_to_data() {
        let code = HammingCode::new(3);
        let data = vec![0, 1, 1, 0];
        let codeword = code.encode(&data);
        let mut corrupted = codeword.clone();
        corrupted[3] ^= 1;
        let recovered = code.decode_to_data(&corrupted);
        assert_eq!(recovered, data);
    }

    #[test]
    fn test_hamming_15_11() {
        let code = HammingCode::new(4);
        assert_eq!(code.n, 15);
        assert_eq!(code.k, 11);
        let data = vec![1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0];
        let codeword = code.encode(&data);
        assert!(code.syndrome(&codeword).iter().all(|&b| b == 0));

        let mut corrupted = codeword.clone();
        corrupted[7] ^= 1;
        let (corrected, pos) = code.decode(&corrupted);
        assert_eq!(corrected, codeword);
        assert_eq!(pos, Some(7));
    }

    #[test]
    fn test_hamming_31_26() {
        let code = HammingCode::new(5);
        assert_eq!(code.n, 31);
        assert_eq!(code.k, 26);
        let data = vec![1; 26];
        let codeword = code.encode(&data);
        let mut corrupted = codeword.clone();
        corrupted[20] ^= 1;
        let (corrected, pos) = code.decode(&corrupted);
        assert_eq!(corrected, codeword);
        assert_eq!(pos, Some(20));
    }

    #[test]
    fn test_hamming_all_zeros() {
        let code = HammingCode::new(3);
        let data = vec![0, 0, 0, 0];
        let codeword = code.encode(&data);
        assert_eq!(codeword, vec![0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_hamming_generator_parity_orthogonality() {
        let code = HammingCode::new(3);
        let g = code.generator_matrix();
        let h = code.parity_check_matrix();
        // G * H^T should be all zeros
        for grow in &g {
            for hrow in &h {
                assert_eq!(gf2::dot(grow, hrow), 0);
            }
        }
    }

    #[test]
    fn test_hamming_all_data_patterns_7_4() {
        let code = HammingCode::new(3);
        for val in 0u8..16 {
            let data: BitVec = (0..4).map(|i| (val >> i) & 1).collect();
            let codeword = code.encode(&data);
            let syndrome = code.syndrome(&codeword);
            assert!(syndrome.iter().all(|&b| b == 0), "Invalid codeword for data {:?}", data);
            let recovered = code.extract_data(&codeword);
            assert_eq!(recovered, data);
        }
    }
}
