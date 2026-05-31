//! Convolutional codes: encoding and Viterbi decoding.

use crate::{Bit, BitVec};

/// A convolutional code defined by generator polynomials.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ConvolutionalCode {
    /// Constraint length K (memory depth + 1).
    pub constraint_length: usize,
    /// Generator polynomials, each as a bit vector (lowest bit = current input).
    pub generators: Vec<BitVec>,
    /// Number of output bits per input bit (rate 1/n).
    pub n_output: usize,
    /// Number of states: 2^(K-1).
    pub num_states: usize,
}

impl ConvolutionalCode {
    /// Create a convolutional code with given constraint length and generator polynomials.
    pub fn new(constraint_length: usize, generators: Vec<BitVec>) -> Self {
        let n_output = generators.len();
        let num_states = 1 << (constraint_length - 1);
        ConvolutionalCode { constraint_length, generators, n_output, num_states }
    }

    /// Get the trellis transition output bits for a given state and input bit.
    /// Returns (new_state, output_bits).
    pub fn transition(&self, state: usize, input: Bit) -> (usize, BitVec) {
        // Shift register: state bits + new input
        let register = (state as u32) | ((input as u32) << (self.constraint_length - 1));
        let mut output = BitVec::with_capacity(self.n_output);
        for gen in &self.generators {
            let mut bit: u8 = 0;
            for (i, &g) in gen.iter().enumerate() {
                if g == 1 {
                    bit ^= ((register >> i) & 1) as u8;
                }
            }
            output.push(bit);
        }
        let new_state = (register >> 1) as usize;
        (new_state, output)
    }

    /// Encode a sequence of input bits.
    pub fn encode(&self, data: &[Bit]) -> BitVec {
        let mut state = 0usize;
        let mut output = BitVec::new();
        for &bit in data {
            let (new_state, bits) = self.transition(state, bit);
            output.extend_from_slice(&bits);
            state = new_state;
        }
        // Terminate: flush with constraint_length - 1 zero bits
        for _ in 0..(self.constraint_length - 1) {
            let (new_state, bits) = self.transition(state, 0);
            output.extend_from_slice(&bits);
            state = new_state;
        }
        output
    }

    /// Viterbi decoding for hard-decision inputs.
    pub fn viterbi_decode(&self, received: &[Bit]) -> BitVec {
        let num_steps = received.len() / self.n_output;
        let num_states = self.num_states;

        // Path metrics
        let mut path_metrics = vec![i64::MAX; num_states];
        path_metrics[0] = 0;

        // Traceback: store the predecessor state and input for each (state, time)
        let mut traceback: Vec<Vec<(usize, Bit)>> = Vec::with_capacity(num_steps);

        for step in 0..num_steps {
            let start = step * self.n_output;
            let end = start + self.n_output;
            if end > received.len() {
                break;
            }
            let recv_chunk = &received[start..end];

            let mut new_metrics = vec![i64::MAX; num_states];
            let mut new_traceback = vec![(0usize, 0u8); num_states];

            for state in 0..num_states {
                if path_metrics[state] == i64::MAX {
                    continue;
                }
                for input in 0u8..2 {
                    let (new_state, expected) = self.transition(state, input);
                    let distance = hamming_distance(recv_chunk, &expected);

                    let new_cost = path_metrics[state] + distance as i64;
                    if new_cost < new_metrics[new_state] {
                        new_metrics[new_state] = new_cost;
                        new_traceback[new_state] = (state, input);
                    }
                }
            }

            path_metrics = new_metrics;
            traceback.push(new_traceback.to_vec());
        }

        // Find the best final state (prefer state 0 for terminated codes)
        let mut best_state = 0;
        let mut best_metric = path_metrics[0];
        for s in 1..num_states {
            if path_metrics[s] < best_metric {
                best_metric = path_metrics[s];
                best_state = s;
            }
        }

        // Traceback
        let mut decoded = BitVec::new();
        let mut state = best_state;
        for step in (0..traceback.len()).rev() {
            let (prev_state, input) = traceback[step][state];
            decoded.push(input);
            state = prev_state;
        }
        decoded.reverse();

        // Remove termination bits
        let data_len = num_steps.saturating_sub(self.constraint_length - 1);
        decoded.truncate(data_len);
        decoded
    }
}

/// Compute Hamming distance between two bit slices.
fn hamming_distance(a: &[Bit], b: &[Bit]) -> usize {
    a.iter().zip(b.iter()).filter(|(&x, &y)| x != y).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_standard_code() -> ConvolutionalCode {
        // Standard rate 1/2, K=3 convolutional code
        // g1 = [1,1,1] (7 octal), g2 = [1,0,1] (5 octal)
        ConvolutionalCode::new(3, vec![
            vec![1, 1, 1],
            vec![1, 0, 1],
        ])
    }

    #[test]
    fn test_convolutional_encode() {
        let code = make_standard_code();
        let data = vec![1, 0, 1, 1];
        let encoded = code.encode(&data);
        // 4 data bits + 2 termination bits = 6 steps, 12 output bits
        assert_eq!(encoded.len(), 12);
    }

    #[test]
    fn test_convolutional_encode_zeros() {
        let code = make_standard_code();
        let data = vec![0, 0, 0];
        let encoded = code.encode(&data);
        // All zeros input with zero termination should produce all zeros
        assert!(encoded.iter().all(|&b| b == 0));
    }

    #[test]
    fn test_viterbi_no_errors() {
        let code = make_standard_code();
        let data = vec![1, 0, 1, 1, 0, 1, 0];
        let encoded = code.encode(&data);
        let decoded = code.viterbi_decode(&encoded);
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_viterbi_single_error() {
        let code = make_standard_code();
        let data = vec![1, 0, 1, 1, 0];
        let mut encoded = code.encode(&data);
        // Introduce a single-bit error
        encoded[3] ^= 1;
        let decoded = code.viterbi_decode(&encoded);
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_viterbi_multiple_errors() {
        let code = make_standard_code();
        let data = vec![1, 1, 0, 1, 0, 0, 1];
        let mut encoded = code.encode(&data);
        // Two errors far enough apart
        encoded[2] ^= 1;
        encoded[8] ^= 1;
        let decoded = code.viterbi_decode(&encoded);
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_viterbi_all_ones() {
        let code = make_standard_code();
        let data = vec![1, 1, 1, 1, 1];
        let encoded = code.encode(&data);
        let decoded = code.viterbi_decode(&encoded);
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_viterbi_long_sequence() {
        let code = make_standard_code();
        let data = vec![1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0];
        let encoded = code.encode(&data);
        let decoded = code.viterbi_decode(&encoded);
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_transition() {
        let code = make_standard_code();
        let (new_state, output) = code.transition(0, 1);
        // Register: 100 (binary), g1=111 -> 1^0^0=1, g2=101 -> 1^0^0=1
        assert_eq!(new_state, 2); // 10 >> 1 = 1... wait, register=4, >>1 = 2
        assert_eq!(output, vec![1, 1]);

        let (new_state, output) = code.transition(2, 0);
        // Register: 010 (state=2, input=0, so register = 2 | 0<<2 = 2 = 010)
        // g1=111: 0^1^0=1, g2=101: 0^0^0=0
        assert_eq!(output, vec![1, 0]);
        assert_eq!(new_state, 1);
    }

    #[test]
    fn test_conv_code_k4() {
        // K=4, rate 1/2
        let code = ConvolutionalCode::new(4, vec![
            vec![1, 1, 1, 1], // 17 octal
            vec![1, 0, 1, 1], // 13 octal
        ]);
        assert_eq!(code.num_states, 8);
        let data = vec![1, 0, 1, 1, 0];
        let encoded = code.encode(&data);
        let decoded = code.viterbi_decode(&encoded);
        assert_eq!(decoded, data);
    }
}
