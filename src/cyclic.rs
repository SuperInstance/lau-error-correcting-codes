//! Cyclic codes: polynomial representation, generator polynomials.

use crate::{Bit, BitVec};

/// A polynomial over GF(2), stored as coefficients from lowest to highest degree.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct BinaryPolynomial {
    /// Coefficients, coeffs[i] is the coefficient of x^i.
    pub coeffs: BitVec,
}

impl BinaryPolynomial {
    /// Create a polynomial from coefficients.
    pub fn new(coeffs: BitVec) -> Self {
        BinaryPolynomial { coeffs }
    }

    /// Create the zero polynomial.
    pub fn zero() -> Self {
        BinaryPolynomial::new(vec![])
    }

    /// Create the constant 1 polynomial.
    pub fn one() -> Self {
        BinaryPolynomial::new(vec![1])
    }

    /// Create x^n (single term).
    pub fn monomial(n: usize) -> Self {
        let mut coeffs = vec![0; n + 1];
        coeffs[n] = 1;
        BinaryPolynomial::new(coeffs)
    }

    /// Degree of the polynomial (None for zero polynomial).
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.iter().rposition(|&c| c == 1)
    }

    /// Evaluate at x=0 (constant term).
    pub fn constant_term(&self) -> Bit {
        self.coeffs.first().copied().unwrap_or(0)
    }

    /// Polynomial addition over GF(2).
    pub fn add(&self, other: &BinaryPolynomial) -> BinaryPolynomial {
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result = vec![0u8; max_len];
        for i in 0..max_len {
            let a = self.coeffs.get(i).copied().unwrap_or(0);
            let b = other.coeffs.get(i).copied().unwrap_or(0);
            result[i] = a ^ b;
        }
        // Trim trailing zeros
        while result.last() == Some(&0) {
            result.pop();
        }
        BinaryPolynomial::new(result)
    }

    /// Polynomial multiplication over GF(2).
    pub fn mul(&self, other: &BinaryPolynomial) -> BinaryPolynomial {
        if self.degree().is_none() || other.degree().is_none() {
            return BinaryPolynomial::zero();
        }
        let deg_a = self.degree().unwrap();
        let deg_b = other.degree().unwrap();
        let mut result = vec![0u8; deg_a + deg_b + 1];
        for i in 0..=deg_a {
            if self.coeffs[i] == 1 {
                for j in 0..=deg_b {
                    result[i + j] ^= other.coeffs[j];
                }
            }
        }
        BinaryPolynomial::new(result)
    }

    /// Polynomial division, returns (quotient, remainder).
    pub fn div_rem(&self, divisor: &BinaryPolynomial) -> (BinaryPolynomial, BinaryPolynomial) {
        let divisor_deg = match divisor.degree() {
            Some(d) => d,
            None => panic!("Division by zero polynomial"),
        };

        let mut remainder = self.coeffs.clone();
        let dividend_deg = match remainder.iter().rposition(|&c| c == 1) {
            Some(d) if d >= divisor_deg => d,
            _ => return (BinaryPolynomial::zero(), BinaryPolynomial::new(remainder)),
        };

        let mut quotient = vec![0u8; dividend_deg - divisor_deg + 1];

        for i in (divisor_deg..=dividend_deg).rev() {
            if remainder.get(i).copied().unwrap_or(0) == 1 {
                quotient[i - divisor_deg] = 1;
                for j in 0..=divisor_deg {
                    let idx = i - divisor_deg + j;
                    if idx < remainder.len() {
                        remainder[idx] ^= divisor.coeffs[j];
                    }
                }
            }
        }

        // Trim remainder
        while remainder.last() == Some(&0) {
            remainder.pop();
        }

        (BinaryPolynomial::new(quotient), BinaryPolynomial::new(remainder))
    }

    /// Check if this polynomial divides other evenly.
    pub fn divides(&self, other: &BinaryPolynomial) -> bool {
        let (_, rem) = other.div_rem(self);
        rem.degree().is_none()
    }

    /// Shift left by n positions (multiply by x^n).
    pub fn shift_left(&self, n: usize) -> BinaryPolynomial {
        let mut coeffs = vec![0u8; n];
        coeffs.extend_from_slice(&self.coeffs);
        BinaryPolynomial::new(coeffs)
    }

    /// Cyclic shift right by 1 for length n.
    pub fn cyclic_shift_right(&self, n: usize) -> BinaryPolynomial {
        let mut coeffs = vec![0u8; self.coeffs.len().max(n)];
        for (i, &c) in self.coeffs.iter().enumerate() {
            if i < n {
                coeffs[(i + n - 1) % n] |= c;
            }
        }
        coeffs.truncate(n);
        BinaryPolynomial::new(coeffs)
    }
}

/// A cyclic code defined by a generator polynomial.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct CyclicCode {
    /// Generator polynomial g(x).
    pub generator: BinaryPolynomial,
    /// Codeword length n.
    pub n: usize,
    /// Data bits k = n - deg(g).
    pub k: usize,
}

impl CyclicCode {
    /// Create a cyclic code with generator g(x) and length n.
    /// g(x) must divide x^n - 1.
    pub fn new(generator: BinaryPolynomial, n: usize) -> Self {
        let k = n - generator.degree().unwrap_or(0);
        CyclicCode { generator, n, k }
    }

    /// Encode data by polynomial multiplication: c(x) = d(x) * g(x).
    pub fn encode(&self, data: &[Bit]) -> BitVec {
        let data_poly = BinaryPolynomial::new(data.to_vec());
        let codeword_poly = data_poly.mul(&self.generator);
        let mut coeffs = codeword_poly.coeffs;
        coeffs.resize(self.n, 0);
        coeffs
    }

    /// Systematic encode: compute remainder of x^(n-k) * d(x) / g(x) and append.
    pub fn encode_systematic(&self, data: &[Bit]) -> BitVec {
        assert_eq!(data.len(), self.k);
        let data_poly = BinaryPolynomial::new(data.to_vec());
        let shifted = data_poly.shift_left(self.n - self.k);
        let (_, remainder) = shifted.div_rem(&self.generator);
        let mut result = data.to_vec();
        // Pad remainder to length n-k
        let mut rem_coeffs = remainder.coeffs.clone();
        rem_coeffs.resize(self.n - self.k, 0);
        result.extend_from_slice(&rem_coeffs);
        result
    }

    /// Compute syndrome: remainder of received(x) / g(x).
    pub fn syndrome(&self, received: &[Bit]) -> BinaryPolynomial {
        let recv_poly = BinaryPolynomial::new(received.to_vec());
        let (_, rem) = recv_poly.div_rem(&self.generator);
        rem
    }

    /// Check if received word is a valid codeword.
    pub fn is_valid(&self, received: &[Bit]) -> bool {
        self.syndrome(received).degree().is_none()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_add() {
        let a = BinaryPolynomial::new(vec![1, 0, 1]); // 1 + x^2
        let b = BinaryPolynomial::new(vec![1, 1]);     // 1 + x
        let sum = a.add(&b);
        assert_eq!(sum, BinaryPolynomial::new(vec![0, 1, 1])); // x + x^2
    }

    #[test]
    fn test_polynomial_mul() {
        let a = BinaryPolynomial::new(vec![1, 1]); // 1 + x
        let b = BinaryPolynomial::new(vec![1, 1]); // 1 + x
        let product = a.mul(&b); // 1 + x^2 (since 2x = 0 in GF(2))
        assert_eq!(product, BinaryPolynomial::new(vec![1, 0, 1]));
    }

    #[test]
    fn test_polynomial_div() {
        // (x^3 + 1) / (x + 1) = x^2 + x + 1 (in GF(2), (x+1)(x^2+x+1) = x^3 + 1)
        let a = BinaryPolynomial::new(vec![1, 0, 0, 1]);
        let b = BinaryPolynomial::new(vec![1, 1]);
        let (q, r) = a.div_rem(&b);
        assert_eq!(q, BinaryPolynomial::new(vec![1, 1, 1]));
        assert_eq!(r, BinaryPolynomial::zero());
    }

    #[test]
    fn test_polynomial_div_with_remainder() {
        let a = BinaryPolynomial::new(vec![1, 1, 0, 1]); // 1 + x + x^3
        let b = BinaryPolynomial::new(vec![1, 1]);        // 1 + x
        let (q, r) = a.div_rem(&b);
        // x^3 + x + 1 = (x+1)(x^2+x) + 1
        assert_eq!(q, BinaryPolynomial::new(vec![0, 1, 1]));
        assert_eq!(r, BinaryPolynomial::new(vec![1]));
    }

    #[test]
    fn test_cyclic_encode() {
        // g(x) = 1 + x + x^3, n = 7 (this is the (7,4) Hamming code as cyclic)
        let g = BinaryPolynomial::new(vec![1, 1, 0, 1]);
        let code = CyclicCode::new(g, 7);
        assert_eq!(code.k, 4);
        let data = vec![1, 0, 1, 0];
        let codeword = code.encode(&data);
        assert!(code.is_valid(&codeword));
    }

    #[test]
    fn test_cyclic_systematic_encode() {
        let g = BinaryPolynomial::new(vec![1, 1, 0, 1]);
        let code = CyclicCode::new(g, 7);
        let data = vec![1, 1, 0, 0];
        let codeword = code.encode_systematic(&data);
        assert_eq!(codeword.len(), 7);
        assert!(code.is_valid(&codeword));
    }

    #[test]
    fn test_cyclic_syndrome() {
        let g = BinaryPolynomial::new(vec![1, 1, 0, 1]);
        let code = CyclicCode::new(g, 7);
        let data = vec![0, 1, 1, 1];
        let codeword = code.encode(&data);
        assert!(code.syndrome(&codeword).degree().is_none());
    }

    #[test]
    fn test_polynomial_degree() {
        assert_eq!(BinaryPolynomial::new(vec![1, 0, 1]).degree(), Some(2));
        assert_eq!(BinaryPolynomial::new(vec![1]).degree(), Some(0));
        assert_eq!(BinaryPolynomial::zero().degree(), None);
    }

    #[test]
    fn test_polynomial_shift() {
        let p = BinaryPolynomial::new(vec![1, 1]); // 1 + x
        let shifted = p.shift_left(2);
        assert_eq!(shifted, BinaryPolynomial::new(vec![0, 0, 1, 1])); // x^2 + x^3
    }
}
