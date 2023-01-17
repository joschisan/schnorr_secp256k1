use super::field;
use super::group;
use super::hash;

const PRIME: [u8; 32] = hex_literal::hex!(
    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
    "FFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"
);
const GROUP_ORDER: [u8; 32] = hex_literal::hex!(
    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE"
    "BAAEDCE6AF48A03BBFD25E8CD0364141"
);

fn hash(a: &[u8; 32], b: &[u8; 32], c: &[u8; 32]) -> [u8; 32] {
    let tag_hash = hex_literal::hex!(
        "7BB52D7A9FEF58323EB1BF7A407DB382"
        "D2F3F2D81BB1224F49FE518F6D48D37C"
    );

    let mut hasher = hash::Sha256::empty();
    hasher.write(&tag_hash);
    hasher.write(&tag_hash);
    hasher.write(a);
    hasher.write(b);
    hasher.write(c);
    hasher.finish()
}

pub fn multiply_by_scalar(point: group::Point, scalar: [u8; 32]) -> group::Point {
    let mut g = group::Point::E;

    let bits = scalar
        .into_iter()
        .flat_map(|x| (0..8).rev().map(move |s| x & (1 << s) != 0));

    for b in bits {
        g = if b { g.double().add(point) } else { g.double() }
    }

    g
}

pub fn public_key(secret_key: [u8; 32]) -> [u8; 32] {
    assert!(secret_key != [0; 32]);
    assert!(secret_key < PRIME);

    multiply_by_scalar(group::GENERATOR, secret_key)
        .affine_x()
        .encode()
}

#[derive(Debug)]
pub enum VerificationError {
    PublicKeyOutOfBounds,
    SignatureOutOfBounds,
    FailedToSolve,
    IsNeutral,
    IsOdd,
    NotEqual,
}

pub fn solve_for_even_y(x: field::Element) -> Result<field::Element, VerificationError> {
    const MAGIC: [u8; 32] = hex_literal::hex!(
        "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
        "FFFFFFFFFFFFFFFFFFFFFFFFBFFFFF0C"
    );

    let c = x.square() * x + 7.into();

    let mut y = field::ONE;

    MAGIC
        .into_iter()
        .flat_map(|x| (0..8).rev().map(move |s| x & (1 << s) != 0))
        .for_each(|b| {
            y = if b { y.square() * c } else { y.square() };
        });

    if y.square() != c {
        return Err(VerificationError::FailedToSolve);
    }

    if !y.is_even() {
        y = y.negative();
    }

    Ok(y)
}

pub fn verify_signature(
    public_key: [u8; 32],
    message: [u8; 32],
    (r, s): ([u8; 32], [u8; 32]),
) -> Result<(), VerificationError> {
    if public_key >= PRIME {
        return Err(VerificationError::PublicKeyOutOfBounds);
    };

    if r >= PRIME || s >= GROUP_ORDER {
        return Err(VerificationError::SignatureOutOfBounds);
    };

    let public_x = field::Element::decode(public_key);
    let public_y = solve_for_even_y(public_x)?;
    let public_point = group::Point::P(public_x, public_y, field::ONE);

    let e = hash(&r, &public_key, &message);

    let g = multiply_by_scalar(group::GENERATOR, s);
    let p = multiply_by_scalar(public_point, e);
    let h = g.add(p.negative());

    if h.is_neutral() {
        return Err(VerificationError::IsNeutral);
    }

    if !h.affine_y().is_even() {
        return Err(VerificationError::IsOdd);
    }

    if h.affine_x().encode() != r {
        return Err(VerificationError::NotEqual);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn solve_for_even_y() {
        let mut g = super::group::GENERATOR;

        for _ in 0..32 {
            let x = g.affine_x();
            let mut y = g.affine_y();

            if !y.is_even() {
                y = y.negative();
            }

            let y_lift = super::solve_for_even_y(x).expect("could not solve for y");

            assert!(y_lift.is_even());
            assert!(y == y_lift);

            g = g.add(super::group::GENERATOR);
        }
    }

    #[test]
    fn public_key() {
        let secret_key: [u8; 32] = hex_literal::hex!(
        "B7E151628AED2A6ABF7158809CF4F3C7"
        "62E7160F38B4DA56A784D9045190CFEF"
        );

        let public_key: [u8; 32] = hex_literal::hex!(
        "DFF1D77F2A671C5F36183726DB2341BE"
        "58FEAE1DA2DECED843240F7B502BA659"
        );

        assert!(super::public_key(secret_key) == public_key);
    }

    #[test]
    fn verify_signature() {
        let public_key: [u8; 32] = hex_literal::hex!(
        "DFF1D77F2A671C5F36183726DB2341BE"
        "58FEAE1DA2DECED843240F7B502BA659"
        );
        let message: [u8; 32] = hex_literal::hex!(
        "243F6A8885A308D313198A2E03707344"
        "A4093822299F31D0082EFA98EC4E6C89"
        );
        let signature_r: [u8; 32] = hex_literal::hex!(
        "6896BD60EEAE296DB48A229FF71DFE07"
        "1BDE413E6D43F917DC8DCF8C78DE3341"
        );
        let signature_s: [u8; 32] = hex_literal::hex!(
        "8906D11AC976ABCCB20B091292BFF4EA"
        "897EFCB639EA871CFA95F6DE339E4B0A"
        );

        let result = super::verify_signature(public_key, message, (signature_r, signature_s));

        assert!(result.is_ok());
    }
}
