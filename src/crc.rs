//! CRC (Cyclic Redundancy Check) computation.

use crate::BitVec;

/// Compute CRC for data using polynomial division.
/// The polynomial is given as a BitVec where the i-th bit is the coefficient of x^i.
/// The polynomial should NOT include the leading x^degree term in the input divisor
/// (it's implicit). However, for full flexibility, we accept the full polynomial.
///
/// Returns the CRC remainder of length = degree(divisor).
pub fn compute_crc(data: &[u8], polynomial: &[u8]) -> BitVec {
    let poly_degree = polynomial.len().saturating_sub(1);
    if poly_degree == 0 {
        return vec![];
    }

    // Append poly_degree zero bits to data
    let mut augmented = data.to_vec();
    augmented.resize(data.len() + poly_degree, 0);

    // Polynomial long division
    for i in 0..data.len() {
        if augmented[i] == 1 {
            for j in 0..polynomial.len() {
                augmented[i + j] ^= polynomial[j];
            }
        }
    }

    // The remainder is the last poly_degree bits
    augmented[data.len()..].to_vec()
}

/// Compute CRC for a byte slice, treating each byte as a sequence of bits (MSB first).
pub fn compute_crc_bytes(data: &[u8], polynomial: &[u8]) -> BitVec {
    let bits = bytes_to_bits(data);
    compute_crc(&bits, polynomial)
}

/// Verify CRC: compute CRC of data+crc and check if remainder is all zeros.
pub fn verify_crc(data_with_crc: &[u8], polynomial: &[u8]) -> bool {
    let remainder = compute_crc(data_with_crc, polynomial);
    remainder.iter().all(|&b| b == 0)
}

/// Common CRC polynomials.
pub mod polynomials {
    use super::BitVec;

    /// CRC-8: x^8 + x^2 + x + 1
    pub fn crc8() -> BitVec {
        vec![1, 0, 0, 0, 0, 0, 1, 1, 1]
    }

    /// CRC-16-CCITT: x^16 + x^12 + x^5 + 1
    pub fn crc16_ccitt() -> BitVec {
        let mut p = vec![0u8; 17];
        p[0] = 1;
        p[5] = 1;
        p[12] = 1;
        p[16] = 1;
        p
    }

    /// CRC-16-IBM: x^16 + x^15 + x^2 + 1
    pub fn crc16_ibm() -> BitVec {
        let mut p = vec![0u8; 17];
        p[0] = 1;
        p[2] = 1;
        p[15] = 1;
        p[16] = 1;
        p
    }

    /// CRC-32: x^32 + x^26 + x^23 + x^22 + x^16 + x^12 + x^11 + x^10 + x^8 + x^7 + x^5 + x^4 + x^2 + x + 1
    pub fn crc32() -> BitVec {
        let mut p = vec![0u8; 33];
        p[0] = 1;
        p[1] = 1;
        p[2] = 1;
        p[4] = 1;
        p[5] = 1;
        p[7] = 1;
        p[8] = 1;
        p[10] = 1;
        p[11] = 1;
        p[12] = 1;
        p[16] = 1;
        p[22] = 1;
        p[23] = 1;
        p[26] = 1;
        p[32] = 1;
        p
    }
}

/// Convert bytes to bits (MSB first).
pub fn bytes_to_bits(bytes: &[u8]) -> BitVec {
    let mut bits = BitVec::with_capacity(bytes.len() * 8);
    for &byte in bytes {
        for i in (0..8).rev() {
            bits.push((byte >> i) & 1);
        }
    }
    bits
}

/// Convert bits to bytes (MSB first), padding with zeros if needed.
pub fn bits_to_bytes(bits: &[u8]) -> Vec<u8> {
    let mut bytes = Vec::new();
    for chunk in bits.chunks(8) {
        let mut byte = 0u8;
        for (i, &bit) in chunk.iter().enumerate() {
            byte |= bit << (7 - i);
        }
        bytes.push(byte);
    }
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crc8_simple() {
        let poly = polynomials::crc8();
        let data = vec![1, 0, 0, 1, 1, 0, 0, 1]; // 0x99 in bits
        let crc = compute_crc(&data, &poly);
        assert_eq!(crc.len(), 8); // CRC-8 produces 8-bit remainder
    }

    #[test]
    fn test_crc_verify() {
        let poly = polynomials::crc8();
        let data = vec![1, 0, 1, 0, 1, 0, 1, 0];
        let crc = compute_crc(&data, &poly);
        let mut full = data.clone();
        full.extend_from_slice(&crc);
        assert!(verify_crc(&full, &poly));
    }

    #[test]
    fn test_crc_detects_error() {
        let poly = polynomials::crc8();
        let data = vec![1, 0, 1, 0, 1, 0, 1, 0];
        let crc = compute_crc(&data, &poly);
        let mut full = data.clone();
        full.extend_from_slice(&crc);
        // Flip a bit
        full[3] ^= 1;
        assert!(!verify_crc(&full, &poly));
    }

    #[test]
    fn test_crc16() {
        let poly = polynomials::crc16_ccitt();
        let data = vec![1, 0, 0, 1, 1, 0, 0, 1];
        let crc = compute_crc(&data, &poly);
        assert_eq!(crc.len(), 16);

        let mut full = data.clone();
        full.extend_from_slice(&crc);
        assert!(verify_crc(&full, &poly));
    }

    #[test]
    fn test_crc32() {
        let poly = polynomials::crc32();
        let data = vec![1; 32];
        let crc = compute_crc(&data, &poly);
        assert_eq!(crc.len(), 32);

        let mut full = data.clone();
        full.extend_from_slice(&crc);
        assert!(verify_crc(&full, &poly));
    }

    #[test]
    fn test_crc_zeros() {
        let poly = polynomials::crc8();
        let data = vec![0; 16];
        let crc = compute_crc(&data, &poly);
        // All zeros should produce all-zero CRC
        assert!(crc.iter().all(|&b| b == 0));
    }

    #[test]
    fn test_bytes_bits_roundtrip() {
        let bytes = vec![0xDE, 0xAD, 0xBE, 0xEF];
        let bits = bytes_to_bits(&bytes);
        assert_eq!(bits.len(), 32);
        let recovered = bits_to_bytes(&bits);
        assert_eq!(recovered, bytes);
    }

    #[test]
    fn test_crc_different_data_different_crc() {
        let poly = polynomials::crc8();
        let data1 = vec![1, 0, 1, 0, 1, 0, 1, 0];
        let data2 = vec![0, 1, 0, 1, 0, 1, 0, 1];
        let crc1 = compute_crc(&data1, &poly);
        let crc2 = compute_crc(&data2, &poly);
        assert_ne!(crc1, crc2);
    }

    #[test]
    fn test_crc_ibm() {
        let poly = polynomials::crc16_ibm();
        let data = vec![1, 0, 1, 1, 0, 0, 1, 1];
        let crc = compute_crc(&data, &poly);
        assert_eq!(crc.len(), 16);

        let mut full = data.clone();
        full.extend_from_slice(&crc);
        assert!(verify_crc(&full, &poly));
    }
}
