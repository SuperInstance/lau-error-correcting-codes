# lau-error-correcting-codes

[![crates.io](https://img.shields.io/badge/crates.io-0.1.0-orange)](https://crates.io/crates/lau-error-correcting-codes)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)
[![docs](https://docs.rs/lau-error-correcting-codes/badge.svg)](https://docs.rs/lau-error-correcting-codes)

**Error-correcting codes in Rust** — parity checks, Hamming codes, linear codes, cyclic codes, Reed-Solomon codes, convolutional codes with Viterbi decoding, CRC, and Shannon's noisy-channel coding theorem.

74 tests · `nalgebra` for linear algebra · `serde` for serialization · GF(2) and GF(2ᵐ) arithmetic

---

## What This Does

When bits travel through a noisy channel, some flip. Error-correcting codes add structured redundancy so the receiver can detect *and fix* errors. This library implements the full classical stack:

- **Parity codes** — single-bit error detection
- **Hamming codes** — single-error correction, double-error detection
- **Linear codes** — generator matrices, parity-check matrices, minimum distance, syndrome decoding
- **Cyclic codes** — polynomial codes over GF(2), generator polynomials
- **Reed-Solomon codes** — burst-error correction over GF(2ᵐ)
- **Convolutional codes** — trellis-based encoding with Viterbi (maximum-likelihood) decoding
- **CRC** — cyclic redundancy checks for error detection
- **Shannon's theorem** — theoretical bounds: channel capacity, entropy, Hamming bound, Singleton bound

Designed for **agent communication reliability** — ensuring error-free message delivery in adversarial or noisy environments.

---

## Key Idea

Claude Shannon proved in 1948 that reliable communication over a noisy channel is possible as long as your transmission rate is below the channel capacity C. This library makes that promise concrete: you pick a code, encode your data, and the decoding algorithm recovers the original message even when the channel corrupts bits.

The hierarchy of codes mirrors a pedagogical journey:

```
Parity → Hamming → Linear → Cyclic → Reed-Solomon → Convolutional
                                          ↕
                                    Shannon's bounds
```

Each level adds capability: parity detects, Hamming corrects, linear codes generalize, cyclic codes add algebraic structure, Reed-Solomon handles bursts, and convolutional codes handle streams.

---

## Install

```toml
[dependencies]
lau-error-correcting-codes = "0.1"
```

Requires **Rust 2021 edition**. Dependencies: `serde`, `nalgebra`.

---

## Quick Start

```rust
use lau_error_correcting_codes::{parity, hamming, linear, crc, shannon};

// 1. Parity check
let data = vec![1, 0, 1, 1, 0, 1, 0];
let parity_bit = parity::even_parity(&data);
assert_eq!(parity_bit, 0); // even number of 1s

// 2. Hamming(7,4) — encode and correct
let hamming = hamming::HammingCode::new(3); // r=3 parity bits → (7,4) code
let encoded = hamming.encode(&vec![1, 0, 1, 1]);
let corrupted = vec![1, 0, 1, 0, 1, 1, 0]; // one bit flipped
let decoded = hamming.decode(&corrupted);
assert_eq!(decoded, vec![1, 0, 1, 1]); // error corrected!

// 3. CRC check
let message = vec![1, 0, 1, 1, 0, 0, 1, 0];
let polynomial = vec![1, 0, 0, 0, 1, 0, 0, 1, 1]; // CRC-8
let checksum = crc::compute_crc(&message, &polynomial);

// 4. Shannon's channel capacity
let capacity = shannon::channel_capacity_bsc(0.01); // BSC with p=0.01
println!("BSC capacity: {:.4} bits/use", capacity);
```

---

## API Reference

### `parity` — Parity Codes

| Function | Description |
|----------|-------------|
| `even_parity(bits)` | Compute even parity bit |
| `odd_parity(bits)` | Compute odd parity bit |
| `check_even_parity(bits)` | Verify even parity of word (including parity bit) |
| `parity_check_matrix(n)` | Construct an (n, n−1) parity-check matrix |

### `hamming` — Hamming Codes

| Type/Method | Description |
|-------------|-------------|
| `HammingCode::new(r)` | Create a Hamming code with `r` parity bits → (2ʳ−1, 2ʳ−r−1) code |
| `encode(data)` | Encode a data word |
| `decode(received)` | Decode via syndrome lookup, correcting up to 1 error |
| `syndrome(received)` | Compute the syndrome vector |
| `generator_matrix()` | Get the generator matrix G |
| `parity_check_matrix()` | Get the parity-check matrix H |

### `linear` — Linear Codes

| Type/Method | Description |
|-------------|-------------|
| `LinearCode::from_generator(generators)` | Define a code by its generator matrix rows |
| `encode(data)` | Linear encoding: c = d·G |
| `decode(received)` | Syndrome decoding using coset leaders |
| `parity_check_matrix()` | Compute H from G |
| `minimum_distance()` | Compute d_min (brute-force for small codes) |
| `code_parameters()` | Return (n, k, d) |

### `cyclic` — Cyclic Codes

| Type/Method | Description |
|-------------|-------------|
| `GF2Polynomial` | Polynomial over GF(2) with arithmetic |
| `CyclicCode::new(n, generator)` | Define a cyclic code by length and generator polynomial |
| `encode(data)` | Systematic encoding via polynomial division |
| `decode(received)` | Meggitt decoder for error correction |

### `reed_solomon` — Reed-Solomon Codes

| Type/Method | Description |
|-------------|-------------|
| `GFElement` | Element of GF(2ᵐ) with full field arithmetic |
| `RSCode::new(n, k, m)` | RS code: n symbols, k data symbols, over GF(2ᵐ) |
| `encode(data)` | Systematic RS encoding |
| `syndrome(received)` | Compute syndromes for error detection |
| `decode(received)` | Berlekamp-Massey + Chien search + Forney correction |

### `convolutional` — Convolutional Codes + Viterbi

| Type/Method | Description |
|-------------|-------------|
| `ConvolutionalCode::new generators, constraint_length)` | Define a code by its generator polynomials |
| `encode(bits)` | Encode a bit stream |
| `viterbi_decode(received)` | Maximum-likelihood decoding via the Viterbi algorithm |

### `crc` — Cyclic Redundancy Check

| Function | Description |
|----------|-------------|
| `compute_crc(data, polynomial)` | Compute CRC checksum via polynomial long division |
| `check_crc(data_with_crc, polynomial)` | Verify CRC (remainder should be zero) |

### `shannon` — Information-Theoretic Bounds

| Function | Description |
|----------|-------------|
| `binary_entropy(p)` | H(p) = −p log₂ p − (1−p) log₂(1−p) |
| `channel_capacity_bsc(p)` | C = 1 − H(p) for a binary symmetric channel |
| `hamming_bound(n, k, t)` | Sphere-packing bound for a (n, k, t) code |
| `singleton_bound(n, k)` | d ≤ n − k + 1 |
| `gilbert_varshamov_bound(n, d)` | Existence bound: codes exist with given n, d |

---

## How It Works

### Encoding Pipeline

1. **Input**: a vector of data bits/symbols
2. **Choose a code**: based on your error model (random → Hamming/linear, burst → Reed-Solomon, streaming → convolutional)
3. **Encode**: add redundant bits via the code's generator matrix or polynomial
4. **Transmit**: send the codeword through the channel
5. **Receive**: get a possibly-corrupted word
6. **Decode**: syndrome-based correction (block codes) or Viterbi traceback (convolutional)

### Decoding Strategies

| Code | Decoding Algorithm | Complexity |
|------|-------------------|------------|
| Hamming | Syndrome lookup table | O(n) |
| Linear | Syndrome + coset leader table | O(2ⁿ⁻ᵏ) precomp, O(n) decode |
| Cyclic | Meggitt decoder | O(n²) |
| Reed-Solomon | Berlekamp-Massey + Chien + Forney | O(n²) |
| Convolutional | Viterbi (dynamic programming on trellis) | O(n · 2ᵏ) |

### GF(2) and GF(2ᵐ) Arithmetic

All binary codes operate over GF(2) — the field {0, 1} with XOR addition and AND multiplication. Reed-Solomon codes operate over extension fields GF(2ᵐ), where elements are polynomials modulo an irreducible polynomial. The library implements full field arithmetic (add, multiply, inverse) for both.

---

## The Math

### Hamming Distance

The **Hamming distance** d(x, y) between two codewords is the number of positions where they differ. A code's **minimum distance** d_min determines its error-correction capability:

> Can correct up to ⌊(d_min − 1)/2⌋ errors
> Can detect up to d_min − 1 errors

### Linear Codes

A linear [n, k, d] code over GF(2) is a k-dimensional subspace of GF(2)ⁿ. It has:

- **Generator matrix** G (k × n): every codeword is c = d·G
- **Parity-check matrix** H ((n−k) × n): every codeword satisfies c·Hᵀ = 0
- **Syndrome** s = r·Hᵀ: s = 0 means no error, s ≠ 0 identifies the error pattern

### Reed-Solomon Codes

An RS(n, k) code over GF(2ᵐ) encodes k data symbols into n codeword symbols using Lagrange interpolation. It achieves the **Singleton bound** with equality: d_min = n − k + 1. This makes it a **Maximum Distance Separable (MDS)** code — the best possible for given n, k.

### Convolutional Codes

Unlike block codes, convolutional codes encode a continuous stream. A rate 1/r code with constraint length K processes one input bit at a time, producing r output bits based on the current bit and K−1 previous bits (the encoder state). The Viterbi algorithm finds the most likely input sequence by finding the shortest path through the **trellis diagram**.

### Shannon's Noisy-Channel Coding Theorem

For a discrete memoryless channel with capacity C:

> For any rate R < C, there exists a code of rate R with arbitrarily low error probability.
> For any rate R > C, the error probability is bounded away from zero.

The capacity of a binary symmetric channel with crossover probability p:

> C = 1 − H(p) = 1 + p log₂ p + (1−p) log₂(1−p)

---

## Test Coverage

| Module | Tests |
|--------|-------|
| `shannon` | 13 |
| `hamming` | 9 |
| `linear` | 9 |
| `cyclic` | 9 |
| `convolutional` | 9 |
| `crc` | 9 |
| `parity` | 9 |
| `reed_solomon` | 7 |
| **Total** | **74** |

---

## License

MIT
