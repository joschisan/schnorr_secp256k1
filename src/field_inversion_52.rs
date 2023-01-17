use super::field;

pub fn invert(x: field::Element) -> field::Element {
    const PRIME: [i64; 5] = [
        0xFFFFEFFFFFC2F,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFF,
        0xFFFFFFFFFFFF,
    ];

    pub const INVERSE_2_POW_728: field::Element = field::Element {
        limbs: [
            0xA9357C4C3E74,
            0x359FCD6E223C,
            0x88010D51154F6,
            0x8F7C917CA2A4C,
            0x45F105208471,
        ],
        magnitude: 0,
    };

    assert!(x != field::ZERO);

    let mut delta = 1;
    let mut f = PRIME.map(|l| l as i128);
    let mut g = x.normalize().limbs.map(|l| l as i128);
    let mut d = [0, 0, 0, 0, 0];
    let mut e = [1, 0, 0, 0, 0];

    let mut transition_matrix;

    for _ in 0..14 {
        (delta, transition_matrix) = update_delta(delta, f[0] as i64, g[0] as i64);

        (f, g) = update_fg(f, g, transition_matrix);
        (d, e) = update_de(d, e, transition_matrix);
    }

    let f = field::Element {
        limbs: [0, 1, 2, 3, 4].map(|i| ((PRIME[i] << 1) + f[i] as i64) as u64),
        magnitude: 2,
    };

    let d = field::Element {
        limbs: [0, 1, 2, 3, 4].map(|i| ((PRIME[i] << 2) + d[i] as i64) as u64),
        magnitude: 3,
    };

    f * d * INVERSE_2_POW_728
}

fn update_delta(mut delta: i64, mut f: i64, mut g: i64) -> (i64, (i64, i64, i64, i64)) {
    //Compute delta and transition matrix t after N divsteps (multiplied by 2^N).

    let (mut u, mut v, mut q, mut r): (i64, i64, i64, i64) = (1, 0, 0, 1); // start with identity matrix

    for _ in 0..52 {
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

    assert!(u.abs() >> 53 == 0);
    assert!(v.abs() >> 53 == 0);
    assert!(q.abs() >> 53 == 0);
    assert!(r.abs() >> 53 == 0);

    (delta, (u, v, q, r))
}

fn update_fg(
    mut f: [i128; 5],
    mut g: [i128; 5],
    transition_matrix: (i64, i64, i64, i64),
) -> ([i128; 5], [i128; 5]) {
    const MASK_52: i128 = 0xFFFFFFFFFFFFF;

    let (u, v, q, r) = transition_matrix;
    let u = u as i128;
    let v = v as i128;
    let q = q as i128;
    let r = r as i128;

    // Start computing t*[f,g].
    let mut cf = u * f[0] + v * g[0];
    let mut cg = q * f[0] + r * g[0];

    // Verify that the bottom 52 bits of the result are zero, and then throw them away.
    assert!(cf & MASK_52 == 0);
    assert!(cg & MASK_52 == 0);

    cf >>= 52;
    cg >>= 52;

    for i in 0..4 {
        // Compute limb i + 1 of t*[f,g], and store it as output limb i shifting down by 52 bits
        cf += u * f[i + 1] + v * g[i + 1];
        cg += q * f[i + 1] + r * g[i + 1];

        f[i] = cf & MASK_52;
        g[i] = cg & MASK_52;

        cf >>= 52;
        cg >>= 52;
    }

    // What remains is limb 5 of t*[f,g]; store it as output limb 4.
    assert!(cf.abs() >> 48 == 0);
    assert!(cg.abs() >> 48 == 0);

    f[4] = cf;
    g[4] = cg;

    (f, g)
}

fn update_de(
    mut d: [i128; 5],
    mut e: [i128; 5],
    transition_matrix: (i64, i64, i64, i64),
) -> ([i128; 5], [i128; 5]) {
    const MASK_48: i128 = 0xFFFFFFFFFFFF;
    const MASK_52: i128 = 0xFFFFFFFFFFFFF;
    // 35 Stellen
    const R: i128 = 0x1000003D1;

    let (u, v, q, r) = transition_matrix;
    let u = u as i128;
    let v = v as i128;
    let q = q as i128;
    let r = r as i128;

    assert!(d[0].abs() >> 53 == 0);
    assert!(d[1].abs() >> 53 == 0);
    assert!(d[2].abs() >> 52 == 0);
    assert!(d[3].abs() >> 52 == 0);
    assert!(d[4].abs() >> 48 == 0);

    assert!(e[0].abs() >> 53 == 0);
    assert!(e[1].abs() >> 53 == 0);
    assert!(e[2].abs() >> 52 == 0);
    assert!(e[3].abs() >> 52 == 0);
    assert!(e[4].abs() >> 48 == 0);

    let mut cd = 0;
    let mut ce = 0;

    for i in 0..4 {
        // Compute limb i
        cd += u * d[i] + v * e[i];
        ce += q * d[i] + r * e[i];

        assert!(cd.abs() >> 108 == 0);
        assert!(ce.abs() >> 108 == 0);

        d[i] = cd & MASK_52;
        e[i] = ce & MASK_52;

        cd >>= 52;
        ce >>= 52;
    }

    cd += u * d[4] + v * e[4];
    ce += q * d[4] + r * e[4];

    assert!(cd.abs() >> 103 == 0);
    assert!(ce.abs() >> 103 == 0);

    d[4] = cd & MASK_48;
    e[4] = ce & MASK_48;

    cd = (cd >> 48) * R;
    ce = (ce >> 48) * R;

    assert!(cd.abs() >> 90 == 0);
    assert!(ce.abs() >> 90 == 0);

    d[0] += cd & MASK_52;
    e[0] += ce & MASK_52;

    cd >>= 52;
    ce >>= 52;

    d[1] += cd;
    e[1] += ce;

    (d, e)
}
