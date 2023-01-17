use crate::field;
use std::u128;

pub fn multiply(a: field::Element, b: field::Element) -> field::Element {
    const MASK_48: u128 = 0xFFFFFFFFFFFF;
    const MASK_52: u128 = 0xFFFFFFFFFFFFF;
    const MASK_64: u128 = 0xFFFFFFFFFFFFFFFF;
    const R: u128 = 0x1000003D10;

    assert!(a.magnitude < 5);
    assert!(b.magnitude < 5);
    a.verify();
    b.verify();

    let a0 = a.limbs[0] as u128;
    let a1 = a.limbs[1] as u128;
    let a2 = a.limbs[2] as u128;
    let a3 = a.limbs[3] as u128;
    let a4 = a.limbs[4] as u128;

    let b0 = b.limbs[0] as u128;
    let b1 = b.limbs[1] as u128;
    let b2 = b.limbs[2] as u128;
    let b3 = b.limbs[3] as u128;
    let b4 = b.limbs[4] as u128;

    //  [... a b c] is a shorthand for ... + a<<104 + b<<52 + c<<0 mod n.
    //  for 0 <= x <= 4, px is a shorthand for sum(a[i]*b[x-i], i=0..x).
    //  for 4 <= x <= 8, px is a shorthand for sum(a[i]*b[x-i], i=(x-4)..4)
    //  Note that [x 0 0 0 0 0] = [x*R].

    let mut d = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0;
    assert!(d >> 114 == 0);
    // [d 0 0 0] = [p3 0 0 0]

    let mut c = a4 * b4;
    assert!(c >> 112 == 0);
    // [c 0 0 0 0 d 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    d += (c & MASK_64) * R;
    c >>= 64;
    assert!(d >> 115 == 0);
    assert!(c >> 48 == 0);
    // [(c<<12) 0 0 0 0 0 d 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    let t3 = d & MASK_52;
    d >>= 52;
    assert!(t3 >> 52 == 0);
    assert!(d >> 63 == 0);
    // [(c<<12) 0 0 0 0 d t3 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    d += a0 * b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4 * b0;
    assert!(d >> 115 == 0);
    // [(c<<12) 0 0 0 0 d t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    d += (c & MASK_64) * (R << 12);
    assert!(d >> 116 == 0);
    // [d t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    let mut t4 = d & MASK_52;
    d >>= 52;
    assert!(t4 >> 52 == 0);
    assert!(d >> 64 == 0);
    // [d t4 t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    let tx = t4 >> 48;
    t4 &= MASK_48;
    assert!(tx >> 4 == 0);
    assert!(t4 >> 48 == 0);
    // [d (t4 + (tx << 48)) t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0] */
    //
    c = a0 * b0;
    assert!(c >> 112 == 0);
    // [d t4+(tx<<48) t3 0 0 c] = [p8 0 0 0 p4 p3 0 0 p0]

    d += a1 * b4 + a2 * b3 + a3 * b2 + a4 * b1;
    assert!(d >> 115 == 0);
    // [d (t4 + (tx << 48)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    let mut u0 = d & MASK_52;
    d >>= 52;
    assert!(u0 >> 52 == 0);
    assert!(d >> 63 == 0);
    // [d u0 (t4 + (tx << 48)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]
    // [d 0 (t4 + (tx << 48) + (u0 << 52)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    u0 = (u0 << 4) | tx;
    assert!(u0 >> 56 == 0);
    // [d 0 t4+(u0<<48) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    c += u0 * (R >> 4);
    assert!(c >> 115 == 0);
    // [d 0 t4 t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    let r0 = c & MASK_52;
    c >>= 52;
    assert!(r0 >> 52 == 0);
    assert!(c >> 61 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 0 p5 p4 p3 0 0 p0]

    c += a0 * b1 + a1 * b0;
    assert!(c >> 114 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 0 p5 p4 p3 0 p1 p0]

    d += a2 * b4 + a3 * b3 + a4 * b2;
    assert!(d >> 114 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    c += (d & MASK_52) * R;
    d >>= 52;
    assert!(c >> 115 == 0);
    assert!(d >> 62 == 0);
    // [d 0 0 t4 t3 0 c r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    let r1 = c & MASK_52;
    c >>= 52;
    assert!(r1 >> 52 == 0);
    assert!(c >> 63 == 0);
    // [d 0 0 t4 t3 c r1 r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    c += a0 * b2 + a1 * b1 + a2 * b0;
    assert!(c >> 114 == 0);
    // [d 0 0 t4 t3 c r1 r0] = [p8 0 p6 p5 p4 p3 p2 p1 p0]

    d += a3 * b4 + a4 * b3;
    assert!(d >> 114 == 0);
    // [d 0 0 t4 t3 c t1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    c += (d & MASK_64) * R;
    d >>= 64;
    assert!(c >> 115 == 0);
    assert!(d >> 50 == 0);
    // [(d<<12) 0 0 0 t4 t3 c r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r2 = c & MASK_52;
    c >>= 52;
    assert!(r2 >> 52 == 0);
    assert!(c >> 63 == 0);
    // [(d<<12) 0 0 0 t4 (t3 + c) r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    c += (d & MASK_64) * (R << 12);
    c += t3;
    assert!(c >> 100 == 0);
    // [t4 c r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r3 = c & MASK_52;
    c >>= 52;
    assert!(r3 >> 52 == 0);
    assert!(c >> 48 == 0);
    // [t4+c r3 r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r4 = (c & MASK_64) + t4;
    assert!(r4 >> 49 == 0);
    // [r4 r3 r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    field::Element {
        limbs: [r0 as u64, r1 as u64, r2 as u64, r3 as u64, r4 as u64],
        magnitude: 0,
    }
}

pub fn square(a: field::Element) -> field::Element {
    const MASK_48: u128 = 0xFFFFFFFFFFFF;
    const MASK_52: u128 = 0xFFFFFFFFFFFFF;
    const MASK_64: u128 = 0xFFFFFFFFFFFFFFFF;
    const R: u128 = 0x1000003D10;

    assert!(a.magnitude < 5);
    a.verify();

    let a0 = a.limbs[0] as u128;
    let a1 = a.limbs[1] as u128;
    let a2 = a.limbs[2] as u128;
    let a3 = a.limbs[3] as u128;
    let a4 = a.limbs[4] as u128;

    //  [... a b c] is a shorthand for ... + a<<104 + b<<52 + c<<0 mod n.
    //  for 0 <= x <= 4, px is a shorthand for sum(a[i]*b[x-i], i=0..x).
    //  for 4 <= x <= 8, px is a shorthand for sum(a[i]*b[x-i], i=(x-4)..4)
    //  Note that [x 0 0 0 0 0] = [x*R].

    let mut d = 2 * (a0 * a3 + a1 * a2);
    assert!(d >> 114 == 0);
    // [d 0 0 0] = [p3 0 0 0]

    let mut c = a4 * a4;
    assert!(c >> 112 == 0);
    // [c 0 0 0 0 d 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    d += (c & MASK_64) * R;
    c >>= 64;
    assert!(d >> 115 == 0);
    assert!(c >> 48 == 0);
    // [(c<<12) 0 0 0 0 0 d 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    let t3 = d & MASK_52;
    d >>= 52;
    assert!(t3 >> 52 == 0);
    assert!(d >> 63 == 0);
    // [(c<<12) 0 0 0 0 d t3 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    d += 2 * (a0 * a4 + a1 * a3) + a2 * a2;
    assert!(d >> 115 == 0);
    // [(c<<12) 0 0 0 0 d t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    d += (c & MASK_64) * (R << 12);
    assert!(d >> 116 == 0);
    // [d t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    let mut t4 = d & MASK_52;
    d >>= 52;
    assert!(t4 >> 52 == 0);
    assert!(d >> 64 == 0);
    // [d t4 t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    let tx = t4 >> 48;
    t4 &= MASK_48;
    assert!(tx >> 4 == 0);
    assert!(t4 >> 48 == 0);
    // [d (t4 + (tx << 48)) t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0] */
    //
    c = a0 * a0;
    assert!(c >> 112 == 0);
    // [d t4+(tx<<48) t3 0 0 c] = [p8 0 0 0 p4 p3 0 0 p0]

    d += 2 * (a1 * a4 + a2 * a3);
    assert!(d >> 115 == 0);
    // [d (t4 + (tx << 48)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    let mut u0 = d & MASK_52;
    d >>= 52;
    assert!(u0 >> 52 == 0);
    assert!(d >> 63 == 0);
    // [d u0 (t4 + (tx << 48)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]
    // [d 0 (t4 + (tx << 48) + (u0 << 52)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    u0 = (u0 << 4) | tx;
    assert!(u0 >> 56 == 0);
    // [d 0 t4+(u0<<48) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    c += u0 * (R >> 4);
    assert!(c >> 115 == 0);
    // [d 0 t4 t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    let r0 = c & MASK_52;
    c >>= 52;
    assert!(r0 >> 52 == 0);
    assert!(c >> 61 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 0 p5 p4 p3 0 0 p0]

    c += 2 * a0 * a1;
    assert!(c >> 114 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 0 p5 p4 p3 0 p1 p0]

    d += 2 * a2 * a4 + a3 * a3;
    assert!(d >> 114 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    c += (d & MASK_52) * R;
    d >>= 52;
    assert!(c >> 115 == 0);
    assert!(d >> 62 == 0);
    // [d 0 0 t4 t3 0 c r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    let r1 = c & MASK_52;
    c >>= 52;
    assert!(r1 >> 52 == 0);
    assert!(c >> 63 == 0);
    // [d 0 0 t4 t3 c r1 r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    c += 2 * a0 * a2 + a1 * a1;
    assert!(c >> 114 == 0);
    // [d 0 0 t4 t3 c r1 r0] = [p8 0 p6 p5 p4 p3 p2 p1 p0]

    d += 2 * a3 * a4;
    assert!(d >> 114 == 0);
    // [d 0 0 t4 t3 c t1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    c += (d & MASK_64) * R;
    d >>= 64;
    assert!(c >> 115 == 0);
    assert!(d >> 50 == 0);
    // [(d<<12) 0 0 0 t4 t3 c r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r2 = c & MASK_52;
    c >>= 52;
    assert!(r2 >> 52 == 0);
    assert!(c >> 63 == 0);
    // [(d<<12) 0 0 0 t4 (t3 + c) r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    c += (d & MASK_64) * (R << 12);
    c += t3;
    assert!(c >> 100 == 0);
    // [t4 c r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r3 = c & MASK_52;
    c >>= 52;
    assert!(r3 >> 52 == 0);
    assert!(c >> 48 == 0);
    // [t4+c r3 r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r4 = (c & MASK_64) + t4;
    assert!(r4 >> 49 == 0);
    // [r4 r3 r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    field::Element {
        limbs: [r0 as u64, r1 as u64, r2 as u64, r3 as u64, r4 as u64],
        magnitude: 0,
    }
}
