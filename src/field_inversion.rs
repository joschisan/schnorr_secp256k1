use super::field;

fn to_i128_offset_62(limbs: [u64; 5]) -> [i128; 5] {
    const MASK_62: u64 = 0x3FFFFFFFFFFFFFFF;

    let l0 = limbs[0] | (limbs[1] << 52);
    let l1 = (limbs[1] >> 10) | (limbs[2] << 42);
    let l2 = (limbs[2] >> 20) | (limbs[3] << 32);
    let l3 = (limbs[3] >> 30) | (limbs[4] << 22);
    let l4 = limbs[4] >> 40;

    [
        l0 & MASK_62,
        l1 & MASK_62,
        l2 & MASK_62,
        l3 & MASK_62,
        l4
    ]
    .map(|l| l as i128)
}

fn to_i64_offset_52(limbs: [i128; 5]) -> [i64; 5] {
    const MASK_52: i64 = 0xFFFFFFFFFFFFF;

    let limbs = limbs.map(|l| l as i64);

    let l0 = limbs[0];
    let l1 = (limbs[1] << 10) + (limbs[0] >> 52);
    let l2 = (limbs[2] << 20) + (limbs[1] >> 42);
    let l3 = (limbs[3] << 30) + (limbs[2] >> 32);
    let l4 = (limbs[4] << 40) + (limbs[3] >> 22);

    [
        l0 & MASK_52,
        l1  & MASK_52,
        l2  & MASK_52,
        l3  & MASK_52,
        l4
    ]
}

pub fn invert(x: field::Element) -> field::Element {
    const PRIME: [u64; 5] = [
        0xFFFFEFFFFFC2F,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFF,
    ];

    const INVERSE_2_POW_728: field::Element = field::Element {
        limbs: [
            0xA9357C4C3E74,
            0x359FCD6E223C,
            0x88010D51154F6,
            0x8F7C917CA2A4C,
            0x45F105208471,
        ],
        magnitude: 0,
    };

    const INVERSE_2_POW_744: field::Element = field::Element {
    limbs: [
        0x223BFB1017899,
        0x54F60359FCD6E,
        0x2A4C88010D511,
        0x84718F7C917CA,
        0xF83445F10520,
    ],
    magnitude: 0x0,
};



    assert!(x != field::ZERO);

    let mut delta = 1;
    let mut f = to_i128_offset_62(PRIME);
    let mut g = to_i128_offset_62(x.normalize().limbs);
    let mut d = [0, 0, 0, 0, 0];
    let mut e = [1, 0, 0, 0, 0];

    let mut transition_matrix;

    for _ in 0..12 {
        (delta, transition_matrix) = update_delta(delta, f[0] as i64, g[0] as i64);

        (f, g) = update_fg(f, g, transition_matrix);
        (d, e) = update_de(d, e, transition_matrix);
    }

    let f = to_i64_offset_52(f);
    let d = to_i64_offset_52(d);

    let f = field::Element {
        limbs: [0, 1, 2, 3, 4].map(|i| ((PRIME[i] << 1) as i64 + f[i]) as u64),
        magnitude: 2,
    };

    let d = field::Element {
        limbs: [0, 1, 2, 3, 4].map(|i| ((PRIME[i] << 2) as i64 + d[i]) as u64),
        magnitude: 3,
    };

    f * d * INVERSE_2_POW_744
}

fn update_delta(mut delta: i64, mut f: i64, mut g: i64) -> (i64, (i64, i64, i64, i64)) {
    //Compute delta and transition matrix t after N divsteps (multiplied by 2^N).

    let (mut u, mut v, mut q, mut r): (i64, i64, i64, i64) = (1, 0, 0, 1); // start with identity matrix

    for _ in 0..62 {
        assert!(f & 1 == 1);
        if delta > 0 && g & 1 == 1 {
            delta = 1 - delta;
            (f, g) = (g, (g - f) >> 1);
            (u, v, q, r) = (q << 1, r << 1, q - u, r - v);
        } else if g & 1 == 1 {
            delta += 1;
            (f, g) = (f, (g + f) >> 1);
            (u, v, q, r) = (u << 1, v << 1, q + u, r + v);
        } else {
            delta += 1;
            (f, g) = (f, g >> 1);
            (u, v, q, r) = (u << 1, v << 1, q, r);
        }
    }

    (delta, (u, v, q, r))
}

fn update_fg(
    mut f: [i128; 5],
    mut g: [i128; 5],
    transition_matrix: (i64, i64, i64, i64),
) -> ([i128; 5], [i128; 5]) {
    const MASK_62: i128 = 0x3FFFFFFFFFFFFFFF;

    let (u, v, q, r) = transition_matrix;
    let u = u as i128;
    let v = v as i128;
    let q = q as i128;
    let r = r as i128;

    // Start computing t*[f,g].
    let mut cf = u * f[0] + v * g[0];
    let mut cg = q * f[0] + r * g[0];

    // Verify that the bottom 62 bits of the result are zero, and then throw them away.
    assert!(cf & MASK_62 == 0);
    assert!(cg & MASK_62 == 0);

    cf >>= 62;
    cg >>= 62;

    for i in 0..4 {
        // Compute limb i + 1 of t*[f,g], and store it as output limb i shifting down by 52 bits
        cf += u * f[i + 1] + v * g[i + 1];
        cg += q * f[i + 1] + r * g[i + 1];

        f[i] = cf & MASK_62;
        g[i] = cg & MASK_62;

        cf >>= 62;
        cg >>= 62;
    }

    // What remains is limb 5 of t*[f,g]; store it as output limb 4.
    assert!(cf.abs() >> 8 == 0);
    assert!(cg.abs() >> 8 == 0);

    f[4] = cf;
    g[4] = cg;

    (f, g)
}

fn update_de(
    mut d: [i128; 5],
    mut e: [i128; 5],
    transition_matrix: (i64, i64, i64, i64),
) -> ([i128; 5], [i128; 5]) {
    const MASK_62: i128 = 0x3FFFFFFFFFFFFFFF;
    const MASK_8: i128 = 0xFF;

    // 35 Stellen
    const R: i128 = 0x1000003D1;

    let (u, v, q, r) = transition_matrix;
    let u = u as i128;
    let v = v as i128;
    let q = q as i128;
    let r = r as i128;

    assert!(d[0].abs() >> 63 == 0);
    assert!(d[1].abs() >> 63 == 0);
    assert!(d[2].abs() >> 62 == 0);
    assert!(d[3].abs() >> 62 == 0);
    assert!(d[4].abs() >> 8 == 0);

    assert!(e[0].abs() >> 63 == 0);
    assert!(e[1].abs() >> 63 == 0);
    assert!(e[2].abs() >> 62 == 0);
    assert!(e[3].abs() >> 62 == 0);
    assert!(e[4].abs() >> 8 == 0);

    let mut cd: i128 = 0;
    let mut ce: i128 = 0;

    for i in 0..4 {
        // Compute limb i
        let xd = u * d[i] + v * e[i];
        let xe = q * d[i] + r * e[i];

        assert!(xd.abs() >> 126 == 0);
        assert!(xe.abs() >> 126 == 0);

        cd += xd;
        ce += xe;

        d[i] = cd & MASK_62;
        e[i] = ce & MASK_62;

        cd >>= 62;
        ce >>= 62;
    }

    cd += u * d[4] + v * e[4];
    ce += q * d[4] + r * e[4];

    assert!(cd.abs() >> 71 == 0);
    assert!(ce.abs() >> 71 == 0);

    d[4] = cd & MASK_8;
    e[4] = ce & MASK_8;

    cd = (cd >> 8) * R;
    ce = (ce >> 8) * R;

    assert!(cd.abs() >> 96 == 0);
    assert!(ce.abs() >> 96 == 0);

    d[0] += cd & MASK_62;
    e[0] += ce & MASK_62;

    cd >>= 62;
    ce >>= 62;

    d[1] += cd;
    e[1] += ce;

    (d, e)
}

#[cfg(test)]
mod tests {
    use crate::field::{self, Element};


    #[test]
    fn conversion() {
        let limbs:[u64;5] = [
            0xA9357C4C3E74,
            0x359FCD6E223C,
            0x88010D51154F6,
            0x8F7C917CA2A4C,
            0x45F105208471
        ];

                println!("{:#X?}",limbs);


        let limbs = super::to_i128_offset_62(limbs);

        let limbs = super::to_i64_offset_52(limbs).map(|l| l as u64);

        println!("{:#X?}",limbs);

       
    }

    #[test]
    fn inverse_744() {

        let mut x = field::ONE;
        for i in 0..744{
            x = x.double();
            x = x.reduce();
        }

        x = x.normalize();


        const INVERSE: [u8; 32] = hex_literal::hex!(
            "FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF"
            "FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2D"
            );

        let mut g = field::ONE;

        let bits = INVERSE
            .into_iter()
            .flat_map(|x| (0..8).rev().map(move |s| x & (1 << s) != 0));

        for b in bits {
            g = if b {
                g.square() * x
            } else {
                g.square()
            }
        }

        g = g.normalize();

        
        println!("{:#X?}", g);

       
    }



}
