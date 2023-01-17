use core::cmp;
use core::convert::From;
use core::ops::{Add, Div, Mul, Sub};

use super::field_inversion;
use super::field_multiplication;

pub const ZERO: Element = Element {
    limbs: [0; 5],
    magnitude: 0,
};

pub const ONE: Element = Element {
    limbs: [1, 0, 0, 0, 0],
    magnitude: 0,
};

#[derive(Debug, Clone, Copy)]
pub struct Element {
    pub limbs: [u64; 5],
    pub magnitude: u64,
}

impl Element {
    pub fn verify(self) {
        assert!((self.limbs[0] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[1] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[2] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[3] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[4] >> (self.magnitude + 49)) == 0);
    }

    pub fn reduce(self) -> Self {
        const MASK_52: u64 = 0xFFFFFFFFFFFFF;
        const MASK_48: u64 = 0xFFFFFFFFFFFF;
        const R: u64 = 0x1000003D1;

        assert!(self.magnitude < 12);
        self.verify();

        let mut limbs = self.limbs;

        limbs[0] += (limbs[4] >> 48) * R;
        limbs[4] &= MASK_48;

        limbs[1] += limbs[0] >> 52;
        limbs[2] += limbs[1] >> 52;
        limbs[3] += limbs[2] >> 52;
        limbs[4] += limbs[3] >> 52;

        limbs[0] &= MASK_52;
        limbs[1] &= MASK_52;
        limbs[2] &= MASK_52;
        limbs[3] &= MASK_52;

        Self {
            limbs,
            magnitude: 0,
        }
    }

    pub fn normalize(self) -> Self {
        const MASK_52: u64 = 0xFFFFFFFFFFFFF;
        const MASK_48: u64 = 0xFFFFFFFFFFFF;
        const R: u64 = 0x1000003D1;
        const P: [u64; 5] = [
            0xFFFFEFFFFFC2F,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFF,
        ];

        // after reduction the value of of a representation is less than 2p
        // since l_4 <= MASK_48 + 2^{12} < 2 * MASK_48 and the bits of l_3
        // that overlap with l_4 are zero
        let mut limbs = self.reduce().limbs;

        // Since there is no overlap in the limbs we can compare the arrays in
        // colexicographical order to determine if the representation is normalized
        let normalized = limbs.iter().rev().lt(P.iter().rev());

        if !normalized {
            // Since a reduced representations value is less then 2p we
            // have to subtract P only once to achieve a normalized representation.
            limbs[0] += R;
            limbs[1] += limbs[0] >> 52;
            limbs[2] += limbs[1] >> 52;
            limbs[3] += limbs[2] >> 52;
            limbs[4] += limbs[3] >> 52;

            // Since the value of a representaion is larger or equal to P
            // and P + R = 2^{256} we carry a 1 to bit 49 of limb 4
            // if it was not set to 1 already. Furthermore the representations
            // value is smaller then 2p + R = 2^{256} + P <= 2^{257}
            assert!(limbs[4] >> 48 == 1);

            limbs[0] &= MASK_52;
            limbs[1] &= MASK_52;
            limbs[2] &= MASK_52;
            limbs[3] &= MASK_52;

            // by zeroing bit 49 we subtract 2^{256} and maintain congruecy modulo P
            limbs[4] &= MASK_48;
        };

        // verify the representation is normalized
        assert!(limbs.iter().rev().lt(P.iter().rev()));

        Self {
            limbs,
            magnitude: 0,
        }
    }

    pub fn decode(bytes: [u8; 32]) -> Self {
        const MASK_52: u64 = 0xFFFFFFFFFFFFF;

        let b3 = u64::from_be_bytes(bytes[..8].try_into().unwrap());
        let b2 = u64::from_be_bytes(bytes[8..16].try_into().unwrap());
        let b1 = u64::from_be_bytes(bytes[16..24].try_into().unwrap());
        let b0 = u64::from_be_bytes(bytes[24..].try_into().unwrap());

        let l0 = b0 & MASK_52;
        let l1 = ((b1 << 12) | b0 >> 52) & MASK_52;
        let l2 = ((b2 << 24) | b1 >> 40) & MASK_52;
        let l3 = ((b3 << 36) | b2 >> 28) & MASK_52;
        let l4 = b3 >> 16;

        Self {
            limbs: [l0, l1, l2, l3, l4],
            magnitude: 0,
        }
    }

    pub fn encode(self) -> [u8; 32] {
        let limbs = self.normalize().limbs;

        let b0 = limbs[0] | (limbs[1] << 52);
        let b1 = (limbs[1] >> 12) | (limbs[2] << 40);
        let b2 = (limbs[2] >> 24) | (limbs[3] << 28);
        let b3 = (limbs[3] >> 36) | (limbs[4] << 16);

        [b3, b2, b1, b0]
            .map(u64::to_be_bytes)
            .concat()
            .try_into()
            .unwrap()
    }

    pub fn is_zero(&self) -> bool {
        self.normalize().limbs == [0, 0, 0, 0, 0]
    }

    pub fn is_even(self) -> bool {
        self.normalize().limbs[0] & 1 == 0
    }

    pub fn double(self) -> Self {
        assert!(self.magnitude < 12);
        self.verify();

        Self {
            limbs: self.limbs.map(|l| l << 1),
            magnitude: self.magnitude + 1,
        }
    }

    pub fn negative(self) -> Self {
        const P: [u64; 5] = [
            0xFFFFEFFFFFC2F,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFF,
        ];

        assert!(self.magnitude < 12);
        self.verify();

        Self {
            limbs: [0, 1, 2, 3, 4].map(|i| (P[i] << (self.magnitude + 1)) - self.limbs[i]),
            magnitude: self.magnitude + 1,
        }
    }

    pub fn square(self) -> Self {
        field_multiplication::square(self)
    }

    pub fn inverse(self) -> Self {
        field_inversion::invert(self)
    }
}

impl PartialEq for Element {
    fn eq(&self, rhs: &Self) -> bool {
        (*self - *rhs).is_zero()
    }
}

impl Eq for Element {}

impl Add for Element {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        assert!(self.magnitude < 12);
        assert!(rhs.magnitude < 12);
        self.verify();
        rhs.verify();

        Self {
            limbs: [0, 1, 2, 3, 4].map(|i| self.limbs[i] + rhs.limbs[i]),
            magnitude: cmp::max(self.magnitude, rhs.magnitude) + 1,
        }
    }
}

impl Sub for Element {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + rhs.negative()
    }
}

impl Mul for Element {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        field_multiplication::multiply(self, rhs)
    }
}

impl Div for Element {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl From<u64> for Element {
    fn from(n: u64) -> Self {
        const MASK_52: u64 = 0xFFFFFFFFFFFFF;
        Self {
            limbs: [n & MASK_52, n >> 52, 0, 0, 0],
            magnitude: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn multiplication() {
        assert!(super::ONE * super::ZERO == super::ZERO);
        assert!(super::ZERO * super::ONE == super::ZERO);
        assert!(super::ONE * super::ONE == super::ONE);
    }

    #[test]
    fn division() {
        const GENERATOR_X: super::Element = super::Element {
            limbs: [
                0x59F2815B16F81,
                0x79829BFCDB2DC,
                0xE28D9055A0629,
                0x5CE870B0779BE,
                0x667EF9DCBBAC,
            ],
            magnitude: 0,
        };

        assert!(super::ONE / super::ONE == super::ONE);
        assert!(GENERATOR_X / GENERATOR_X == super::ONE);
    }
}
