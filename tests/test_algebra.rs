use isogeny::quaternion::algebra::{QuatAlg, QuatElem};
use num_bigint::BigInt;
use isogeny::quaternion::algebra::BigIntAlg;
use isogeny::quaternion::bigint_crypto::CryptoInt;


type TestAlg = QuatAlg<BigInt>;
type TestElem = QuatElem<BigInt>;

const LIMBS: usize = 4;
type CryptoTestAlg = QuatAlg<CryptoInt<LIMBS>>;
type CryptoTestElem = QuatElem<CryptoInt<LIMBS>>;


#[test]
fn test_quat_elem_initialization_and_conj() {
    let elem = TestElem::new_i32(2, 1, -3, 4, -5);

    assert_eq!(elem.denom, BigInt::from(2));
    assert_eq!(elem.coord[0], BigInt::from(1));
    assert_eq!(elem.coord[1], BigInt::from(-3));

    let conj = elem.conj();

    assert_eq!(conj.coord[0], elem.coord[0]);
    assert_eq!(conj.coord[1], -elem.coord[1].clone());
    assert_eq!(conj.coord[2], -elem.coord[2].clone());
    assert_eq!(conj.coord[3], -elem.coord[3].clone());
}

#[test]
fn test_quat_elem_addition() {
    let a = TestElem::new_i32(3, 1, 2, 3, 4);
    let b = TestElem::new_i32(2, 5, -1, 0, 2);

    let sum = a.add(&b);

    // a = (1 + 2i + 3j + 4ij) / 3
    // b = (5 -  i      + 2ij) / 2
    // a + b = (17 + i + 6j + 14ij) / 6

    assert_eq!(sum.denom, BigInt::from(6));
    assert_eq!(sum.coord[0], BigInt::from(17));
    assert_eq!(sum.coord[1], BigInt::from(1));
    assert_eq!(sum.coord[2], BigInt::from(6));
    assert_eq!(sum.coord[3], BigInt::from(14));
}

#[test]
fn test_quat_elem_multiplication() {
    // p = 11
    let alg = TestAlg::new(BigInt::from(11));

    let a = TestElem::new_i32(2, 1, 2, 3, 4); // (1 + 2i + 3j + 4ij) / 2
    let b = TestElem::new_i32(3, 5, 6, 7, 8); // (5 + 6i + 7j + 8ij) / 3

    let prod = a.mul(&b, &alg);
    // -295/3 - 14/3*i + 5*j + 4*k
    // but not normalized!

    assert_eq!(prod.denom, BigInt::from(6));
    assert_eq!(prod.coord[0], BigInt::from(-590));
    assert_eq!(prod.coord[1], BigInt::from(-28));
    assert_eq!(prod.coord[2], BigInt::from(30));
    assert_eq!(prod.coord[3], BigInt::from(24));
}

#[test]
fn test_quat_elem_norm() {
    let alg = TestAlg::new(BigInt::from(11));
    let a = TestElem::new_i32(2, 1, -2, 3, -4);
    // a = (1 - 2i + 3j - 4ij)/2; p = 11
    // nrd(a) = 280/4

    let (num, denom) = a.norm(&alg);

    assert_eq!(num, BigInt::from(70));
    assert_eq!(denom, BigInt::from(1));
}

#[test]
fn test_quat_elem_normalize() {
    let mut a = TestElem::new_i32(-10, -4, 6, -8, 2);
    a.normalize();

    // (-4 + 6i - 8j + 2ij)/(-10) = (2 - 3i + 4j - ij)/5

    assert_eq!(a.denom, BigInt::from(5));
    assert_eq!(a.coord[0], BigInt::from(2));
    assert_eq!(a.coord[1], BigInt::from(-3));
    assert_eq!(a.coord[2], BigInt::from(4));
    assert_eq!(a.coord[3], BigInt::from(-1));
}

// Same tests but with bigint_crypto backend

#[test]
fn test_crypto_quat_elem_initialization_and_conj() {
    let elem = CryptoTestElem::new_i32(2, 1, -3, 4, -5);

    assert_eq!(elem.denom, CryptoInt::from_i32(2));
    assert_eq!(elem.coord[0], CryptoInt::from_i32(1));
    assert_eq!(elem.coord[1], CryptoInt::from_i32(-3));

    let conj = elem.conj();

    assert_eq!(conj.coord[0], elem.coord[0]);
    assert_eq!(conj.coord[1], -elem.coord[1].clone());
    assert_eq!(conj.coord[2], -elem.coord[2].clone());
    assert_eq!(conj.coord[3], -elem.coord[3].clone());
}

#[test]
fn test_crypto_quat_elem_addition() {
    let a = CryptoTestElem::new_i32(3, 1, 2, 3, 4);
    let b = CryptoTestElem::new_i32(2, 5, -1, 0, 2);

    let sum = a.add(&b);

    assert_eq!(sum.denom, CryptoInt::from_i32(6));
    assert_eq!(sum.coord[0], CryptoInt::from_i32(17));
    assert_eq!(sum.coord[1], CryptoInt::from_i32(1));
    assert_eq!(sum.coord[2], CryptoInt::from_i32(6));
    assert_eq!(sum.coord[3], CryptoInt::from_i32(14));
}

#[test]
fn test_crypto_quat_elem_multiplication() {
    // p = 11
    let alg = CryptoTestAlg::new(CryptoInt::from_i32(11));

    let a = CryptoTestElem::new_i32(2, 1, 2, 3, 4);
    let b = CryptoTestElem::new_i32(3, 5, 6, 7, 8);

    let prod = a.mul(&b, &alg);

    assert_eq!(prod.denom, CryptoInt::from_i32(6));
    assert_eq!(prod.coord[0], CryptoInt::from_i32(-590));
    assert_eq!(prod.coord[1], CryptoInt::from_i32(-28));
    assert_eq!(prod.coord[2], CryptoInt::from_i32(30));
    assert_eq!(prod.coord[3], CryptoInt::from_i32(24));
}

#[test]
fn test_crypto_quat_elem_norm() {
    let alg = CryptoTestAlg::new(CryptoInt::from_i32(11));
    let a = CryptoTestElem::new_i32(2, 1, 1, 1, 1);

    let (num, denom) = a.norm(&alg);

    assert_eq!(num, CryptoInt::from_i32(6));
    assert_eq!(denom, CryptoInt::from_i32(1));
}

#[test]
fn test_crypto_quat_elem_normalize() {
    let mut a = CryptoTestElem::new_i32(-10, -4, 6, -8, 2);

    a.normalize();

    assert_eq!(a.denom, CryptoInt::from_i32(5));
    assert_eq!(a.coord[0], CryptoInt::from_i32(2));
    assert_eq!(a.coord[1], CryptoInt::from_i32(-3));
    assert_eq!(a.coord[2], CryptoInt::from_i32(4));
    assert_eq!(a.coord[3], CryptoInt::from_i32(-1));
}
