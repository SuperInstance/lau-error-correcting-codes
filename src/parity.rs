//! Parity codes: even/odd parity, parity check matrices.

use crate::{Bit, BitVec, gf2};

/// Compute even parity bit for a slice of bits.
/// Returns 1 if the number of 1s is odd (to make total even), 0 otherwise.
pub fn even_parity(data: &[Bit]) -> Bit {
    data.iter().fold(0u8, |acc, &b| acc ^ b)
}

/// Compute odd parity bit for a slice of bits.
/// Returns 1 if the number of 1s is even (to make total odd), 0 otherwise.
pub fn odd_parity(data: &[Bit]) -> Bit {
    1 - even_parity(data)
}

/// Append an even parity bit to the data.
pub fn encode_even_parity(data: &[Bit]) -> BitVec {
    let mut encoded = data.to_vec();
    encoded.push(even_parity(data));
    encoded
}

/// Append an odd parity bit to the data.
pub fn encode_odd_parity(data: &[Bit]) -> BitVec {
    let mut encoded = data.to_vec();
    encoded.push(odd_parity(data));
    encoded
}

/// Check if even parity is valid.
pub fn check_even_parity(data_with_parity: &[Bit]) -> bool {
    if data_with_parity.is_empty() {
        return true;
    }
    let (data, parity) = data_with_parity.split_at(data_with_parity.len() - 1);
    even_parity(data) == parity[0]
}

/// Check if odd parity is valid.
pub fn check_odd_parity(data_with_parity: &[Bit]) -> bool {
    if data_with_parity.is_empty() {
        return true;
    }
    let (data, parity) = data_with_parity.split_at(data_with_parity.len() - 1);
    odd_parity(data) == parity[0]
}

/// A parity check matrix H such that H * x^T = 0 for valid codewords.
/// Stored as a vector of rows, where each row is a bit vector.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ParityCheckMatrix {
    /// Rows of the parity check matrix.
    pub rows: Vec<BitVec>,
    /// Number of columns (codeword length).
    pub n: usize,
}

impl ParityCheckMatrix {
    /// Create a new parity check matrix from rows.
    pub fn new(rows: Vec<BitVec>) -> Self {
        let n = rows.first().map(|r| r.len()).unwrap_or(0);
        ParityCheckMatrix { rows, n }
    }

    /// Compute syndrome: H * received^T.
    pub fn syndrome(&self, received: &[Bit]) -> BitVec {
        self.rows.iter().map(|row| gf2::dot(row, received)).collect()
    }

    /// Check if a received word is a valid codeword (syndrome is all zeros).
    pub fn is_valid(&self, received: &[Bit]) -> bool {
        self.syndrome(received).iter().all(|&b| b == 0)
    }
}

/// Build a simple single-parity-check matrix for n-bit codewords (1 parity bit).
pub fn single_parity_check_matrix(n: usize) -> ParityCheckMatrix {
    let row = vec![1u8; n];
    ParityCheckMatrix::new(vec![row])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_even_parity_all_zeros() {
        assert_eq!(even_parity(&[0, 0, 0, 0]), 0);
    }

    #[test]
    fn test_even_parity_odd_ones() {
        assert_eq!(even_parity(&[1, 0, 1, 0, 1]), 1);
    }

    #[test]
    fn test_even_parity_even_ones() {
        assert_eq!(even_parity(&[1, 0, 1, 0]), 0);
    }

    #[test]
    fn test_odd_parity_all_zeros() {
        assert_eq!(odd_parity(&[0, 0, 0, 0]), 1);
    }

    #[test]
    fn test_encode_check_even_parity() {
        let data = vec![1, 0, 1, 1];
        let encoded = encode_even_parity(&data);
        assert_eq!(encoded.len(), 5);
        assert!(check_even_parity(&encoded));
        // Flip a bit
        let mut corrupted = encoded.clone();
        corrupted[2] ^= 1;
        assert!(!check_even_parity(&corrupted));
    }

    #[test]
    fn test_encode_check_odd_parity() {
        let data = vec![1, 0, 1, 1];
        let encoded = encode_odd_parity(&data);
        assert!(check_odd_parity(&encoded));
    }

    #[test]
    fn test_parity_check_matrix_valid() {
        let data = vec![1, 0, 1, 0];
        let encoded = encode_even_parity(&data);
        let h = single_parity_check_matrix(5);
        assert!(h.is_valid(&encoded));
    }

    #[test]
    fn test_parity_check_matrix_invalid() {
        let data = vec![1, 0, 1, 0];
        let mut encoded = encode_even_parity(&data);
        encoded[1] ^= 1;
        let h = single_parity_check_matrix(5);
        assert!(!h.is_valid(&encoded));
        let syndrome = h.syndrome(&encoded);
        assert_eq!(syndrome, vec![1]);
    }

    #[test]
    fn test_custom_parity_check_matrix() {
        // H = [[1,1,0],[0,1,1]]
        let h = ParityCheckMatrix::new(vec![
            vec![1, 1, 0],
            vec![0, 1, 1],
        ]);
        // Valid codeword: [0,0,0] -> syndrome [0,0]
        assert!(h.is_valid(&[0, 0, 0]));
        // Valid codeword: [1,1,1] -> syndrome [0,0]
        assert!(h.is_valid(&[1, 1, 1]));
        // Invalid: [1,0,0] -> syndrome [1,0]
        assert!(!h.is_valid(&[1, 0, 0]));
    }
}
