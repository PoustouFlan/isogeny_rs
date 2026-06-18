use isogeny::quaternion::algebra::{BigIntAlg, QuatAlg, QuatElem};
use isogeny::quaternion::lattice::QuatLattice;
use isogeny::quaternion::matrix::MatrixUtils;
use num_bigint::BigInt;


#[inline]
fn b_zero() -> BigInt { <BigInt as BigIntAlg>::zero() }
#[inline]
fn b(val: i32) -> BigInt { <BigInt as BigIntAlg>::from_i32(val) }

#[test]
fn test_lattice_equal() {
    let mut lat = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    lat.basis = MatrixUtils::identity();
    cmp.basis = MatrixUtils::identity();
    assert!(QuatLattice::equal(&lat, &cmp));

    lat.denom = b(5);
    cmp.denom = b(4);
    assert!(!QuatLattice::equal(&lat, &cmp));

    lat.denom = b(1);
    cmp.denom = b(-1);
    assert!(QuatLattice::equal(&lat, &cmp));

    lat.denom = b(3);
    cmp.denom = b(3);
    assert!(QuatLattice::equal(&lat, &cmp));

    lat.basis[0][0] = b(1);
    lat.basis[0][3] = b(-1);
    lat.basis[1][1] = b(-2);
    lat.basis[2][2] = b(1);
    lat.basis[2][1] = b(1);
    lat.basis[3][3] = b(-3);
    lat.denom = b(6);
    lat.hnf();

    cmp.basis = lat.basis.clone();
    cmp.denom = b(6);
    assert!(QuatLattice::equal(&lat, &cmp));

    cmp.denom = b(-7);
    assert!(!QuatLattice::equal(&lat, &cmp));

    cmp.denom = b(6);
    cmp.basis[3][3] = b(165);
    assert!(!QuatLattice::equal(&lat, &cmp));
}

#[test]
fn test_lattice_inclusion() {
    let mut lat = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    lat.basis = MatrixUtils::identity();
    cmp.basis = MatrixUtils::identity();
    assert!(QuatLattice::inclusion(&lat, &cmp));

    lat.denom = b(5);
    cmp.denom = b(4);
    assert!(!QuatLattice::inclusion(&lat, &cmp));

    lat.denom = b(1);
    cmp.denom = b(3);
    assert!(QuatLattice::inclusion(&lat, &cmp));

    lat.denom = b(3);
    cmp.denom = b(3);
    assert!(QuatLattice::inclusion(&lat, &cmp));

    lat.basis[0][0] = b(1);
    lat.basis[0][3] = b(-1);
    lat.basis[1][1] = b(-2);
    lat.basis[2][2] = b(1);
    lat.basis[2][1] = b(1);
    lat.basis[3][3] = b(-3);
    lat.denom = b(6);
    lat.hnf();

    cmp.basis = lat.basis.clone();
    cmp.denom = b(6);
    assert!(QuatLattice::inclusion(&lat, &cmp));

    cmp.denom = b(12);
    assert!(QuatLattice::inclusion(&lat, &cmp));

    cmp.denom = b(6);
    cmp.basis[3][3] = b(165);
    assert!(!QuatLattice::inclusion(&lat, &cmp));
}

#[test]
fn test_lattice_reduce_denom() {
    let mut lat = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    let s = 15;
    for i in 0..4 {
        for j in 0..4 {
            lat.basis[i][j] = b((i as i32 + j as i32) * s);
            cmp.basis[i][j] = b(i as i32 + j as i32);
        }
    }
    lat.denom = b(4 * s);
    cmp.denom = b(4);

    let mut red = lat.clone();
    red.reduce_denom();

    assert!(MatrixUtils::equal(&red.basis, &cmp.basis));
    assert_eq!(red.denom, cmp.denom);

    lat.reduce_denom();
    assert!(MatrixUtils::equal(&lat.basis, &cmp.basis));
    assert_eq!(lat.denom, cmp.denom);
}

#[test]
fn test_lattice_conjugate_without_hnf() {
    let mut lat = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    lat.basis = MatrixUtils::zero();
    lat.basis[0][0] = b(4);
    lat.basis[0][3] = b(1);
    lat.basis[1][1] = b(-2);
    lat.basis[2][2] = b(-1);
    lat.basis[2][1] = b(-1);
    lat.basis[3][3] = b(-3);
    lat.denom = b(6);

    cmp.basis = MatrixUtils::zero();
    cmp.basis[0][0] = b(4);
    cmp.basis[0][3] = b(1);
    cmp.basis[1][1] = b(2);
    cmp.basis[2][2] = b(1);
    cmp.basis[2][1] = b(1);
    cmp.basis[3][3] = b(3);
    cmp.denom = b(6);

    lat.hnf();

    let mut conj = lat.conjugate_without_hnf();
    conj.hnf();
    cmp.hnf();

    assert!(QuatLattice::equal(&conj, &cmp));

    let mut conj_conj = conj.conjugate_without_hnf();
    conj_conj.hnf();
    assert!(QuatLattice::equal(&conj_conj, &lat));
}

#[test]
fn test_lattice_dual_without_hnf() {
    let mut lat = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    lat.basis = MatrixUtils::zero();
    lat.basis[0][0] = b(1);
    lat.basis[0][3] = b(-1);
    lat.basis[1][1] = b(-2);
    lat.basis[2][2] = b(1);
    lat.basis[2][1] = b(1);
    lat.basis[3][3] = b(-3);
    lat.denom = b(6);

    cmp.basis = MatrixUtils::zero();
    cmp.basis[0][0] = b(6);
    cmp.basis[1][1] = b(3);
    cmp.basis[2][2] = b(6);
    cmp.basis[3][3] = b(2);
    cmp.denom = b(1);

    lat.hnf();

    let mut dual = lat.dual_without_hnf();
    dual.hnf();
    cmp.hnf();

    assert!(QuatLattice::equal(&dual, &cmp));
    assert!(!QuatLattice::equal(&dual, &lat));

    let mut dual_dual = dual.dual_without_hnf();
    dual_dual.hnf();
    assert!(QuatLattice::equal(&dual_dual, &lat));
}

#[test]
fn test_lattice_add() {
    let mut lat1 = QuatLattice::<BigInt>::zero();
    let mut lat2 = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    lat1.basis = MatrixUtils::zero();
    lat2.basis = MatrixUtils::zero();
    cmp.basis = MatrixUtils::zero();

    lat1.basis[0][0] = b(44);
    lat1.basis[0][2] = b(3);
    lat1.basis[0][3] = b(32);
    lat2.basis[0][0] = b(1);
    cmp.basis[0][0] = b(2);
    cmp.basis[0][2] = b(1);
    lat1.basis[1][1] = b(5);
    lat2.basis[1][1] = b(2);
    cmp.basis[1][1] = b(1);
    lat1.basis[2][2] = b(3);
    lat2.basis[2][2] = b(1);
    cmp.basis[2][2] = b(1);
    lat1.basis[3][3] = b(1);
    lat2.basis[3][3] = b(3);
    cmp.basis[3][3] = b(3);

    lat1.denom = b(4);
    lat2.denom = b(6);
    cmp.denom = b(12);

    let sum = QuatLattice::add(&lat1, &lat2);
    assert!(MatrixUtils::equal(&sum.basis, &cmp.basis));
    assert_eq!(sum.denom, cmp.denom);

    lat1.basis = MatrixUtils::zero();
    lat2.basis = MatrixUtils::zero();

    lat1.basis[0][0] = b(4);
    lat1.basis[0][2] = b(3);
    lat2.basis[0][0] = b(1);
    lat2.basis[0][3] = b(-1);
    lat1.basis[1][1] = b(5);
    lat2.basis[1][1] = b(-2);
    lat1.basis[2][2] = b(3);
    lat2.basis[2][2] = b(1);
    lat2.basis[2][1] = b(1);
    lat1.basis[3][3] = b(7);
    lat2.basis[3][3] = b(-3);
    lat1.denom = b(4);
    lat2.denom = b(6);

    let sum2 = QuatLattice::add(&lat1, &lat2);
    assert!(MatrixUtils::equal(&sum2.basis, &cmp.basis));
    assert_eq!(sum2.denom, cmp.denom);

    cmp.basis = lat2.basis.clone();
    cmp.denom = lat2.denom.clone();
    cmp.hnf();

    let sum_self = QuatLattice::add(&lat2, &lat2);
    assert!(MatrixUtils::equal(&sum_self.basis, &cmp.basis));
    assert_eq!(sum_self.denom, cmp.denom);
}

#[test]
fn test_lattice_intersect() {
    let mut lat1 = QuatLattice::<BigInt>::zero();
    let mut lat2 = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    lat1.basis = MatrixUtils::zero();
    lat2.basis = MatrixUtils::zero();
    cmp.basis = MatrixUtils::zero();

    lat1.basis[0][0] = b(4);
    lat1.basis[0][2] = b(3);
    lat2.basis[0][0] = b(1);
    lat2.basis[0][3] = b(-1);
    lat1.basis[1][1] = b(5);
    lat2.basis[1][1] = b(-2);
    lat1.basis[2][2] = b(3);
    lat2.basis[2][2] = b(1);
    lat2.basis[2][1] = b(1);
    lat1.basis[3][3] = b(7);
    lat2.basis[3][3] = b(-3);
    lat1.denom = b(4);
    lat2.denom = b(6);
    lat1.hnf();
    lat2.hnf();

    cmp.basis[0][0] = b(2);
    cmp.basis[0][2] = b(1);
    cmp.basis[1][1] = b(10);
    cmp.basis[2][2] = b(3);
    cmp.basis[3][3] = b(7);
    cmp.denom = b(2);

    let inter = QuatLattice::intersect(&lat1, &lat2);
    assert!(QuatLattice::equal(&inter, &cmp));

    let inter2 = QuatLattice::intersect(&lat2, &lat1);
    assert!(QuatLattice::equal(&inter2, &cmp));

    cmp.basis = lat1.basis.clone();
    cmp.denom = lat1.denom.clone();
    let inter_self = QuatLattice::intersect(&lat1, &lat1);
    assert!(QuatLattice::equal(&inter_self, &cmp));
}

#[test]
fn test_lattice_alg_elem_mul() {
    let mut lat = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();
    let alg = QuatAlg::new(b(23));

    let elem = QuatElem::new_i32(2, 3, 4, -1, 0);

    lat.basis = MatrixUtils::zero();
    lat.basis[0][0] = b(11);
    lat.basis[1][1] = b(-13);
    lat.basis[2][2] = b(15);
    lat.basis[3][3] = b(-4);
    lat.basis[0][1] = b(2);
    lat.basis[1][3] = b(-1);
    lat.denom = b(5);
    lat.hnf();

    let prod = QuatLattice::alg_elem_mul(&lat, &elem, &alg);

    cmp.basis = MatrixUtils::zero();
    cmp.basis[0][0] = b(33);
    cmp.basis[1][0] = b(44);
    cmp.basis[2][0] = b(-11);
    cmp.basis[3][0] = b_zero();
    cmp.basis[0][1] = b(27 - 4 * 13);
    cmp.basis[1][1] = b(36 + 3 * 13);
    cmp.basis[2][1] = b(-9);
    cmp.basis[3][1] = b(-13);
    cmp.basis[0][2] = b(15 * 23);
    cmp.basis[1][2] = b_zero();
    cmp.basis[2][2] = b(45);
    cmp.basis[3][2] = b(-60);
    cmp.basis[0][3] = b(-4);
    cmp.basis[1][3] = b(3 + 23 * 4);
    cmp.basis[2][3] = b(16);
    cmp.basis[3][3] = b(-1 + 4 * 3);
    cmp.denom = b(10);
    cmp.hnf();

    assert!(QuatLattice::equal(&cmp, &prod));

    let prod2 = QuatLattice::alg_elem_mul(&lat, &elem, &alg);
    assert!(QuatLattice::equal(&cmp, &prod2));
}

#[test]
fn test_lattice_mul() {
    let mut lat1 = QuatLattice::<BigInt>::zero();
    let mut lat2 = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();
    let alg = QuatAlg::new(b(19));

    lat1.basis = MatrixUtils::zero();
    lat2.basis = MatrixUtils::zero();
    cmp.basis = MatrixUtils::zero();

    lat1.basis[0][0] = b(44);
    lat1.basis[0][2] = b(3);
    lat1.basis[0][3] = b(32);
    lat2.basis[0][0] = b(1);
    cmp.basis[0][0] = b(1);
    lat1.basis[1][1] = b(5);
    lat2.basis[1][1] = b(2);
    cmp.basis[1][1] = b(1);
    lat1.basis[2][2] = b(3);
    lat2.basis[2][2] = b(1);
    cmp.basis[2][2] = b(1);
    lat1.basis[3][3] = b(1);
    lat2.basis[3][3] = b(3);
    cmp.basis[3][3] = b(1);
    lat1.denom = b(4);
    lat2.denom = b(6);
    cmp.denom = b(24);

    let prod = QuatLattice::mul(&lat1, &lat2, &alg);
    assert!(MatrixUtils::equal(&prod.basis, &cmp.basis));
    assert_eq!(prod.denom, cmp.denom);

    lat1.basis = MatrixUtils::zero();
    lat2.basis = MatrixUtils::zero();

    lat1.basis[0][0] = b(4);
    lat1.basis[0][2] = b(3);
    lat2.basis[0][0] = b(1);
    lat2.basis[0][3] = b(-1);
    lat1.basis[1][1] = b(5);
    lat2.basis[1][1] = b(-2);
    lat1.basis[2][2] = b(3);
    lat2.basis[2][2] = b(1);
    lat2.basis[2][1] = b(1);
    lat1.basis[3][3] = b(7);
    lat2.basis[3][3] = b(-3);
    lat1.denom = b(4);
    lat2.denom = b(6);

    let prod2 = QuatLattice::mul(&lat1, &lat2, &alg);
    assert!(MatrixUtils::equal(&prod2.basis, &cmp.basis));
    assert_eq!(prod2.denom, cmp.denom);

    cmp.basis = MatrixUtils::zero();
    cmp.basis[0][0] = b(1);
    cmp.basis[1][1] = b(1);
    cmp.basis[2][2] = b(1);
    cmp.basis[3][3] = b(1);
    cmp.denom = b(36);

    let prod_self = QuatLattice::mul(&lat2, &lat2, &alg);
    assert!(MatrixUtils::equal(&prod_self.basis, &cmp.basis));
    assert_eq!(prod_self.denom, cmp.denom);
}

#[test]
fn test_lattice_contains() {
    // let alg = QuatAlg::new(b(103));
    let mut lat = QuatLattice::<BigInt>::zero();

    lat.basis = MatrixUtils::zero();
    lat.basis[0][0] = b(4);
    lat.basis[0][2] = b(3);
    lat.basis[1][1] = b(5);
    lat.basis[2][2] = b(3);
    lat.basis[3][3] = b(7);
    lat.denom = b(4);

    let x = QuatElem::new_i32(3, 1, -2, 26, 9);
    assert!(lat.contains(&x).is_none());

    lat.basis = MatrixUtils::zero();
    lat.basis[0][0] = b(1);
    lat.basis[0][3] = b(-1);
    lat.basis[1][1] = b(-2);
    lat.basis[2][2] = b(1);
    lat.basis[2][1] = b(1);
    lat.basis[3][3] = b(-3);
    lat.denom = b(6);
    lat.hnf();

    let coord_res = lat.contains(&x);
    assert!(coord_res.is_some());
    let coord = coord_res.unwrap();
    assert_eq!(coord[0], b(2));
    assert_eq!(coord[1], b(-2));
    assert_eq!(coord[2], b(52));
    assert_eq!(coord[3], b(6));
}

#[test]
fn test_lattice_index() {
    let mut sublat = QuatLattice::<BigInt>::zero();
    let mut overlat = QuatLattice::<BigInt>::zero();

    sublat.basis = MatrixUtils::zero();
    overlat.basis = MatrixUtils::identity();
    overlat.denom = b(2);

    sublat.basis[0][0] = b(2);
    sublat.basis[0][1] = b_zero();
    sublat.basis[0][2] = b(1);
    sublat.basis[0][3] = b_zero();

    sublat.basis[1][0] = b_zero();
    sublat.basis[1][1] = b(4);
    sublat.basis[1][2] = b(2);
    sublat.basis[1][3] = b(3);

    sublat.basis[2][0] = b_zero();
    sublat.basis[2][1] = b_zero();
    sublat.basis[2][2] = b(1);
    sublat.basis[2][3] = b_zero();

    sublat.basis[3][0] = b_zero();
    sublat.basis[3][1] = b_zero();
    sublat.basis[3][2] = b_zero();
    sublat.basis[3][3] = b(1);

    sublat.denom = b(2);

    let index = QuatLattice::index(&sublat, &overlat);
    assert_eq!(index, b(8));
}

#[test]
fn test_lattice_hnf() {
    let mut lat = QuatLattice::<BigInt>::zero();
    let mut cmp = QuatLattice::<BigInt>::zero();

    lat.basis = MatrixUtils::zero();
    cmp.basis = MatrixUtils::zero();

    lat.basis[0][0] = b(1);
    lat.basis[0][3] = b(-1);
    lat.basis[1][1] = b(-2);
    lat.basis[2][2] = b(1);
    lat.basis[2][1] = b(1);
    lat.basis[3][3] = b(-3);

    cmp.basis[0][0] = b(1);
    cmp.basis[1][1] = b(2);
    cmp.basis[2][2] = b(1);
    cmp.basis[3][3] = b(3);

    cmp.denom = b(6);
    lat.denom = b(6);

    lat.hnf();

    assert!(MatrixUtils::equal(&lat.basis, &cmp.basis));
    assert_eq!(lat.denom, cmp.denom);
}

#[test]
fn test_lattice_gram() {
    let mut lattice = QuatLattice::<BigInt>::zero();
    let alg = QuatAlg::new(b(103));

    lattice.basis = MatrixUtils::zero();
    lattice.basis[0][0] = b(202);
    lattice.basis[1][1] = b(202);
    lattice.basis[2][2] = b(1);
    lattice.basis[3][3] = b(1);
    lattice.basis[0][2] = b(158);
    lattice.basis[0][3] = b(53);
    lattice.basis[1][2] = b(149);
    lattice.basis[1][3] = b(158);
    lattice.denom = b(2);

    let gram = lattice.gram(&alg);

    let elem1 = QuatElem::new_i32(2, 360, 149, 1, 0);
    let elem2 = QuatElem::new_i32(2, 53, 360, 0, 1);

    let vec1_opt = lattice.contains(&elem1);
    let vec2_opt = lattice.contains(&elem2);
    assert!(vec1_opt.is_some());
    assert!(vec2_opt.is_some());

    let vec1 = vec1_opt.unwrap();
    let vec2 = vec2_opt.unwrap();

    let mut elem2_conj = elem2.clone();
    elem2_conj = elem2_conj.conj();
    let prod = elem1.mul(&elem2_conj, &alg);

    let norm1 = (prod.coord[0].clone() * b(2)) / prod.denom.clone();

    let mut vec1_gram = [b_zero(), b_zero(), b_zero(), b_zero()];
    for i in 0..4 {
        let mut sum = b_zero();
        for j in 0..4 {
            sum = sum + gram[i][j].clone() * vec1[j].clone();
        }
        vec1_gram[i] = sum;
    }

    let mut norm2 = b_zero();
    for i in 0..4 {
        norm2 = norm2 + vec1_gram[i].clone() * vec2[i].clone();
    }

    norm2 = norm2 / lattice.denom.clone();
    norm2 = norm2 / lattice.denom.clone();

    assert_eq!(norm1, norm2);
}
