const K: [u32; 64] = [
    0x428A2F98, 0x71374491, 0xB5C0FBCF, 0xE9B5DBA5, 0x3956C25B, 0x59F111F1, 0x923F82A4, 0xAB1C5ED5,
    0xD807AA98, 0x12835B01, 0x243185BE, 0x550C7DC3, 0x72BE5D74, 0x80DEB1FE, 0x9BDC06A7, 0xC19BF174,
    0xE49B69C1, 0xEFBE4786, 0xFC19DC6, 0x240CA1CC, 0x2DE92C6F, 0x4A7484AA, 0x5CB0A9DC, 0x76F988DA,
    0x983E5152, 0xA831C66D, 0xB00327C8, 0xBF597FC7, 0xC6E00BF3, 0xD5A79147, 0x6CA6351, 0x14292967,
    0x27B70A85, 0x2E1B2138, 0x4D2C6DFC, 0x53380D13, 0x650A7354, 0x766A0ABB, 0x81C2C92E, 0x92722C85,
    0xA2BFE8A1, 0xA81A664B, 0xC24B8B70, 0xC76C51A3, 0xD192E819, 0xD6990624, 0xF40E3585, 0x106AA070,
    0x19A4C116, 0x1E376C08, 0x2748774C, 0x34B0BCB5, 0x391C0CB3, 0x4ED8AA4A, 0x5B9CCA4F, 0x682E6FF3,
    0x748F82EE, 0x78A5636F, 0x84C87814, 0x8CC70208, 0x90BEFFFA, 0xA4506CEB, 0xBEF9A3F7, 0xC67178F2,
];

const H: [u32; 8] = [
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A, 0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
];

pub struct Sha256 {
    state: [u32; 8],
    buffer: [u8; 64],
    n_buffer: usize,
    n_rounds: u64,
}

impl Sha256 {
    pub fn empty() -> Self {
        Self {
            state: H,
            buffer: [0u8; 64],
            n_buffer: 0,
            n_rounds: 0,
        }
    }

    pub fn write(&mut self, mut bytes: &[u8]) {
        if self.n_buffer + bytes.len() >= 64 {
            let (b1, b2) = bytes.split_at(64 - self.n_buffer);
            bytes = b2;
            self.buffer[self.n_buffer..].copy_from_slice(b1);
            update_state(&mut self.state, &self.buffer);
            self.n_buffer = 0;
            self.n_rounds += 1;
        }

        while bytes.len() >= 64 {
            let (b1, b2) = bytes.split_at(64);
            bytes = b2;
            update_state(&mut self.state, b1.try_into().unwrap());
            self.n_rounds += 1;
        }

        assert!(self.n_buffer + bytes.len() <= 64);

        self.buffer[self.n_buffer..self.n_buffer + bytes.len()].copy_from_slice(bytes);
        self.n_buffer += bytes.len();
    }

    pub fn finish(mut self) -> [u8; 32] {
        let n_written_bits = self.n_rounds * 512 + (self.n_buffer as u64) * 8;
        let n_padding = if self.n_buffer < 56 { 55 } else { 119 } - self.n_buffer;
        let padding = &[0; 63][..n_padding];

        self.write(&[128]);
        self.write(padding);
        self.write(&n_written_bits.to_be_bytes());

        assert!(self.n_buffer == 0);

        unsafe { std::mem::transmute::<[u32; 8], [u8; 32]>(self.state.map(u32::to_be)) }
    }
}

fn update_state(state: &mut [u32; 8], bytes: &[u8; 64]) {
    let mut w: [u32; 64] = [0; 64];
    let byte_blocks = unsafe { std::mem::transmute::<&[u8; 64], &[[u8; 4]; 16]>(bytes) };
    let words: [u32; 16] = byte_blocks.map(u32::from_be_bytes);
    w[0..16].copy_from_slice(words.as_slice());

    for i in 16..64 {
        let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
        let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16]
            .wrapping_add(s0)
            .wrapping_add(w[i - 7])
            .wrapping_add(s1);
    }

    let mut h = *state;
    for i in 0..64 {
        let ch = (h[4] & h[5]) ^ (!h[4] & h[6]);
        let ma = (h[0] & h[1]) ^ (h[0] & h[2]) ^ (h[1] & h[2]);
        let s0 = h[0].rotate_right(2) ^ h[0].rotate_right(13) ^ h[0].rotate_right(22);
        let s1 = h[4].rotate_right(6) ^ h[4].rotate_right(11) ^ h[4].rotate_right(25);
        let t0 = h[7]
            .wrapping_add(s1)
            .wrapping_add(ch)
            .wrapping_add(K[i])
            .wrapping_add(w[i]);
        let t1 = s0.wrapping_add(ma);

        h[7] = h[6];
        h[6] = h[5];
        h[5] = h[4];
        h[4] = h[3].wrapping_add(t0);
        h[3] = h[2];
        h[2] = h[1];
        h[1] = h[0];
        h[0] = t0.wrapping_add(t1);
    }

    for i in 0..8 {
        state[i] = state[i].wrapping_add(h[i]);
    }
}

#[cfg(test)]
mod tests {
    fn digest(data: &[u8]) -> [u8; 32] {
        let mut sha256 = super::Sha256::empty();
        sha256.write(data);
        sha256.finish()
    }

    /// Run empty test input from FIPS 180-2
    #[test]
    fn sha256_nist_empty() {
        let digest = digest(&[]);
        let hash: [u8; 32] = [
            0xE3, 0xB0, 0xC4, 0x42, 0x98, 0xFC, 0x1C, 0x14, 0x9A, 0xFB, 0xF4, 0xC8, 0x99, 0x6F,
            0xB9, 0x24, 0x27, 0xAE, 0x41, 0xE4, 0x64, 0x9B, 0x93, 0x4C, 0xA4, 0x95, 0x99, 0x1B,
            0x78, 0x52, 0xB8, 0x55,
        ];

        assert!(digest == hash);
    }

    /// Run abc test from FIPS 180-2
    #[test]
    fn sha256_nist_abc() {
        let digest = digest("abc".as_bytes());
        let hash: [u8; 32] = [
            0xBA, 0x78, 0x16, 0xBF, 0x8F, 0x1, 0xCF, 0xEA, 0x41, 0x41, 0x40, 0xDE, 0x5D, 0xAE,
            0x22, 0x23, 0xB0, 0x3, 0x61, 0xA3, 0x96, 0x17, 0x7A, 0x9C, 0xB4, 0x10, 0xFF, 0x61,
            0xF2, 0x0, 0x15, 0xAD,
        ];

        assert!(digest == hash);
    }

    /// Run two-block test from FIPS 180-2
    #[test]
    fn sha256_nist_two_blocks() {
        let digest = digest("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq".as_bytes());
        let hash: [u8; 32] = [
            0x24, 0x8D, 0x6A, 0x61, 0xD2, 0x6, 0x38, 0xB8, 0xE5, 0xC0, 0x26, 0x93, 0xC, 0x3E, 0x60,
            0x39, 0xA3, 0x3C, 0xE4, 0x59, 0x64, 0xFF, 0x21, 0x67, 0xF6, 0xEC, 0xED, 0xD4, 0x19,
            0xDB, 0x6, 0xC1,
        ];

        assert!(digest == hash);
    }

    /// Run large input test (1,000,000 x a) from FIPS 180-2
    #[test]
    fn sha256_nist_large_input() {
        let input_str = std::iter::repeat("a").take(1_000_000).collect::<String>();
        let digest = digest(input_str.as_bytes());
        let hash: [u8; 32] = [
            0xCD, 0xC7, 0x6E, 0x5C, 0x99, 0x14, 0xFB, 0x92, 0x81, 0xA1, 0xC7, 0xE2, 0x84, 0xD7,
            0x3E, 0x67, 0xF1, 0x80, 0x9A, 0x48, 0xA4, 0x97, 0x20, 0xE, 0x4, 0x6D, 0x39, 0xCC, 0xC7,
            0x11, 0x2C, 0xD0,
        ];
        assert!(digest == hash);
    }
}
