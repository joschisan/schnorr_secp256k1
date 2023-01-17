use super::field;

const THREE_HALFS: field::Element = field::Element {
    limbs: [
        0xFFFFF7FFFFE19,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFF,
        0x7FFFFFFFFFFF,
    ],
    magnitude: 0,
};

const GENERATOR_X: field::Element = field::Element {
    limbs: [
        0x2815B16F81798,
        0xDB2DCE28D959F,
        0xE870B07029BFC,
        0xBBAC55A06295C,
        0x79BE667EF9DC,
    ],
    magnitude: 0,
};

const GENERATOR_Y: field::Element = field::Element {
    limbs: [
        0x7D08FFB10D4B8,
        0x48A68554199C4,
        0xE1108A8FD17B4,
        0xC4655DA4FBFC0,
        0x483ADA7726A3,
    ],
    magnitude: 0,
};

pub const GENERATOR: Point = Point::P(GENERATOR_X, GENERATOR_Y, field::ONE);

#[derive(Debug, Clone, Copy)]
pub enum Point {
    E,
    P(field::Element, field::Element, field::Element),
}

impl Point {
    fn verify(self) {
        if let Point::P(x, y, z) = self {
            let z2 = z.square();
            let x = x / z2;
            let y = y / (z2 * z);
            assert!(y.square() == x * x * x + 7.into());
        }
    }

    pub fn affine_x(self) -> field::Element {
        match self {
            Self::E => panic!("affine_x() has been called on the neutral element"),
            Self::P(x, y, z) => x / z.square(),
        }
    }

    pub fn affine_y(self) -> field::Element {
        match self {
            Self::E => panic!("affine_x() has been called on the neutral element"),
            Self::P(x, y, z) => y / (z.square() * z),
        }
    }

    pub fn is_neutral(self) -> bool {
        match self {
            Point::E => true,
            Point::P(..) => false,
        }
    }

    pub fn negative(self) -> Self {
        match self {
            Self::E => Self::E,
            Self::P(x, y, z) => Self::P(x, y.negative(), z),
        }
    }

    pub fn double(self) -> Self {
        match self {
            Point::E => Point::E,
            Point::P(x, y, z) => {
                if y.is_zero() {
                    return Point::E;
                }

                /* Formula used:
                 * L = (3/2) * X^2
                 * S = Y^2
                 * T = -X*S
                 * RX = L^2 + 2*T
                 * RY = -(L*(RX + T) + S^2)
                 * RZ = Y*Z
                 */

                let l = THREE_HALFS * x.square();
                let s = y.square();
                let t = (x * s).negative();
                let rx = l.square() + t.double();
                let ry = (l * (rx + t) + s.square()).negative();
                let rz = y * z;

                Point::P(rx, ry, rz)
            }
        }
    }

    pub fn add(self, rhs: Self) -> Self {
        match (self, rhs) {
            (Point::E, Point::E) => Point::E,
            (Point::E, Point::P(..)) => rhs,
            (Point::P(..), Point::E) => self,
            (Point::P(ax, ay, az), Point::P(bx, by, bz)) => {
                let az2 = az.square();
                let bz2 = bz.square();

                let ax = ax * bz2;
                let ay = ay * bz2 * bz;
                let bx = bx * az2;
                let by = by * az2 * az;

                let h = bx - ax;
                let i = ay - by;

                if h.is_zero() {
                    if i.is_zero() {
                        return self.double();
                    } else {
                        assert!(ay == by.negative());
                        return Self::E;
                    }
                }

                let h2 = h.square().negative();
                let h3 = h2 * h;
                let t = ax * h2;

                let rx = i.square() + h3 + t.double();
                let ry = i * (t + rx) + h3 * ay;
                let rz = az * bz * h;

                Self::P(rx, ry, rz)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn verify_generator() {
        super::GENERATOR.verify();
    }

    #[test]
    fn verify_double() {
        let mut g = super::GENERATOR;
        for _ in 0..32 {
            g = g.double();
            g.verify();
        }
    }

    #[test]
    fn verify_add() {
        let mut g = super::GENERATOR;
        for _ in 0..32 {
            g = g.add(super::GENERATOR);
            g.verify();
        }
    }

    #[test]
    fn little_fermat() {
        const CURVE_ORDER: [u8; 32] = hex_literal::hex!(
        "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE"
        "BAAEDCE6AF48A03BBFD25E8CD0364141"
        );

        let mut g = super::Point::E;

        let bits = CURVE_ORDER
            .into_iter()
            .flat_map(|x| (0..8).rev().map(move |s| x & (1 << s) != 0));

        for b in bits {
            g = if b {
                g.double().add(super::GENERATOR)
            } else {
                g.double()
            }
        }

        assert!(g.is_neutral());
    }
}
