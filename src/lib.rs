//! # lau-error-correcting-codes
//!
//! Error-correcting codes for detecting and correcting errors in data transmission.
//! Designed for agent communication reliability — ensuring error-free message delivery.
//!
//! ## Modules
//! - `parity` — Even/odd parity, parity check matrices
//! - `hamming` — Hamming codes with syndrome decoding
//! - `cyclic` — Cyclic codes with polynomial representation
//! - `reed_solomon` — Reed-Solomon basics
//! - `linear` — Linear codes, generator/parity matrices, minimum distance
//! - `convolutional` — Convolutional codes with Viterbi decoding
//! - `crc` — Cyclic redundancy check
//! - `shannon` — Shannon's noisy channel coding theorem bounds

pub mod parity;
pub mod hamming;
pub mod cyclic;
pub mod reed_solomon;
pub mod linear;
pub mod convolutional;
pub mod crc;
pub mod shannon;

/// Represents a bit as a u8 (0 or 1).
pub type Bit = u8;

/// A vector of bits.
pub type BitVec = Vec<u8>;

/// Galois Field GF(2) arithmetic helpers.
pub mod gf2 {
    use crate::Bit;

    /// Add two bits in GF(2) (XOR).
    #[inline]
    pub fn add(a: Bit, b: Bit) -> Bit {
        a ^ b
    }

    /// Multiply two bits in GF(2) (AND).
    #[inline]
    pub fn mul(a: Bit, b: Bit) -> Bit {
        a & b
    }

    /// Dot product of two bit vectors in GF(2).
    pub fn dot(a: &[Bit], b: &[Bit]) -> Bit {
        a.iter().zip(b.iter()).fold(0u8, |acc, (&x, &y)| acc ^ (x & y))
    }

    /// Add two bit vectors in GF(2).
    pub fn vec_add(a: &[Bit], b: &[Bit]) -> Vec<Bit> {
        a.iter().zip(b.iter()).map(|(&x, &y)| x ^ y).collect()
    }

    /// Multiply a bit vector by a scalar in GF(2).
    pub fn vec_scale(v: &[Bit], s: Bit) -> Vec<Bit> {
        v.iter().map(|&x| x & s).collect()
    }
}

/// Galois Field GF(2^m) element for Reed-Solomon codes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct GFElement {
    /// The element value (polynomial coefficient representation).
    pub val: u16,
    /// The field order m (GF(2^m)).
    pub m: u8,
}

impl GFElement {
    pub fn new(val: u16, m: u8) -> Self {
        GFElement { val, m }
    }

    /// Irreducible polynomial for GF(2^m). Returns the reduction polynomial.
    fn reducing_poly(m: u8) -> u16 {
        match m {
            3 => 0b1011,     // x^3 + x + 1
            4 => 0b10011,    // x^4 + x + 1
            8 => 0b100011011, // x^8 + x^4 + x^3 + x^2 + 1 (AES)
            _ => 0b10011,    // default x^4 + x + 1
        }
    }

    /// The field size: 2^m.
    pub fn field_size(&self) -> u16 {
        1u16 << self.m
    }

    /// Add two GF(2^m) elements.
    pub fn add(&self, other: &GFElement) -> GFElement {
        GFElement::new(self.val ^ other.val, self.m)
    }

    /// Subtract (same as add in GF(2^m)).
    pub fn sub(&self, other: &GFElement) -> GFElement {
        self.add(other)
    }

    /// Multiply two GF(2^m) elements.
    pub fn mul(&self, other: &GFElement) -> GFElement {
        let m = self.m;
        let reduce = Self::reducing_poly(m);
        let mut result: u16 = 0;
        let mut a = self.val;
        let b = other.val;
        for i in 0..m {
            if (b >> i) & 1 == 1 {
                result ^= a;
            }
            a <<= 1;
            if a & (1 << m) != 0 {
                a ^= reduce;
            }
        }
        GFElement::new(result, m)
    }

    /// Compute multiplicative inverse in GF(2^m) via exhaustive search.
    pub fn inv(&self) -> Option<GFElement> {
        if self.val == 0 {
            return None;
        }
        let one = GFElement::new(1, self.m);
        for v in 1..self.field_size() {
            let candidate = GFElement::new(v, self.m);
            if self.mul(&candidate) == one {
                return Some(candidate);
            }
        }
        None
    }

    /// Divide two GF(2^m) elements.
    pub fn div(&self, other: &GFElement) -> Option<GFElement> {
        other.inv().map(|inv| self.mul(&inv))
    }

    /// Power: self^exp.
    pub fn pow(&self, mut exp: u32) -> GFElement {
        let m = self.m;
        let mut result = GFElement::new(1, m);
        let mut base = *self;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exp >>= 1;
        }
        result
    }
}

impl std::fmt::Display for GFElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GF(2^{})[{}]", self.m, self.val)
    }
}
