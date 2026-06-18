// tests/test_hnf.rs

use isogeny::quaternion::algebra::BigIntAlg;
use isogeny::quaternion::hnf::*;
use isogeny::quaternion::matrix::{Mat4x4, MatrixUtils};
use num_bigint::BigInt;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

// Helper functions to unambiguously target BigIntAlg methods
#[inline]
fn b_zero() -> BigInt { <BigInt as BigIntAlg>::zero() }
#[inline]
fn b_one() -> BigInt { <BigInt as BigIntAlg>::one() }

/// Helper function to strictly verify if a matrix is in SQISign's Row Hermite Normal Form
fn assert_is_hnf(mat: &Mat4x4<BigInt>) {
    let zero = b_zero();
    
    for i in 0..4 {
        // 1. Must be upper triangular (lower-left is 0)
        for j in 0..i {
            assert_eq!(
                mat[i][j], zero,
                "HNF Failed: Not upper triangular at [{}][{}]. Value: {:?}", i, j, mat[i][j]
            );
        }

        // 2. Pivot (diagonal) must be strictly positive
        assert!(
            mat[i][i] > zero,
            "HNF Failed: Diagonal pivot not positive at [{}][{}]. Value: {:?}", i, i, mat[i][i]
        );

        // 3. Elements TO THE RIGHT of the pivot must be >= 0 and strictly < ROW pivot
        for j in (i + 1)..4 {
            assert!(
                mat[i][j] >= zero,
                "HNF Failed: Element to right is negative at [{}][{}]. Value: {:?}", i, j, mat[i][j]
            );
            assert!(
                mat[i][j] < mat[i][i],
                "HNF Failed: Element to right >= row pivot at [{}][{}]. Value: {:?} >= {:?}",
                i, j, mat[i][j], mat[i][i]
            );
        }
    }
}

#[test]
fn test_mod_not_zero() {
    let modulo = BigInt::from_i32(3);
    
    assert_eq!(mod_not_zero(&BigInt::from_i32(5), &modulo), BigInt::from_i32(2));
    assert_eq!(mod_not_zero(&BigInt::from_i32(6), &modulo), BigInt::from_i32(3));
    assert_eq!(mod_not_zero(&BigInt::from_i32(-1), &modulo), BigInt::from_i32(2));
}

#[test]
fn test_centered_mod() {
    let modulo = BigInt::from_i32(3);
    
    assert_eq!(centered_mod(&BigInt::from_i32(5), &modulo), BigInt::from_i32(-1));
    assert_eq!(centered_mod(&BigInt::from_i32(4), &modulo), BigInt::from_i32(1));
    assert_eq!(centered_mod(&BigInt::from_i32(6), &modulo), b_zero());
}

#[test]
fn test_xgcd_with_u_not_0() {
    let a = b_zero();
    let b = BigInt::from_i32(5);
    
    let (d, u, v) = xgcd_with_u_not_0(&a, &b);
    
    assert_eq!(d, BigInt::from_i32(5));
    assert!(!BigIntAlg::is_zero(&u), "xgcd_with_u_not_0 failed to guarantee u != 0");
    
    let res = a.clone() * u.clone() + b.clone() * v.clone();
    assert_eq!(res, d, "Bezout identity failed");
}

#[test]
fn test_mat_4x4_inv_with_det_as_denom() {
    
    let mut mat = MatrixUtils::zero::<BigInt>();
    mat[0][0] = BigInt::from_i32(1);
    mat[1][1] = BigInt::from_i32(4);
    mat[2][2] = BigInt::from_i32(1);
    mat[3][3] = BigInt::from_i32(1);
    mat[0][1] = BigInt::from_i32(2);
    mat[1][0] = BigInt::from_i32(3);

    let mut adj = MatrixUtils::zero::<BigInt>();
    let det = mat_4x4_inv_with_det_as_denom(Some(&mut adj), &mat);

    assert_eq!(det, BigInt::from_i32(-2));

    let mut prod = MatrixUtils::zero::<BigInt>();
    for i in 0..4 {
        for j in 0..4 {
            let mut sum = b_zero();
            for k in 0..4 {
                sum = sum + adj[i][k].clone() * mat[k][j].clone();
            }
            prod[i][j] = sum;
        }
    }

    let mut expected = MatrixUtils::zero::<BigInt>();
    for i in 0..4 {
        expected[i][i] = det.clone();
    }

    assert!(MatrixUtils::equal(&prod, &expected), "Adjugate multiplication property failed");
}

#[test]
fn test_mat_4xn_hnf_mod_core_exact_vectors() {

    let mut generators = vec![[b_zero(), b_zero(), b_zero(), b_zero()]; 8];
    generators[0][0] = BigInt::from(4);
    generators[2][0] = BigInt::from(3);
    generators[4][0] = BigInt::from(1);
    generators[7][0] = BigInt::from(-1);
    generators[1][1] = BigInt::from(5);
    generators[5][1] = BigInt::from(-2);
    generators[2][2] = BigInt::from(3);
    generators[6][2] = BigInt::from(1);
    generators[5][2] = BigInt::from(1);
    generators[3][3] = BigInt::from(7);
    generators[7][3] = BigInt::from(-3);

    let det = BigInt::from(1);
    let mut hnf = MatrixUtils::zero::<BigInt>();
    
    mat_4xn_hnf_mod_core(&mut hnf, &generators, &det);
    
    let identity = MatrixUtils::identity::<BigInt>();
    assert!(MatrixUtils::equal(&identity, &hnf), "Identity vector test failed");

    let mut generators = vec![[b_zero(), b_zero(), b_zero(), b_zero()]; 8];
    generators[4][0] = BigInt::from(438);
    generators[4][1] = BigInt::from(400);
    generators[4][2] = BigInt::from(156);
    generators[4][3] = BigInt::from(-2);
    generators[5][0] = BigInt::from(-400);
    generators[5][1] = BigInt::from(438);
    generators[5][2] = BigInt::from(2);
    generators[5][3] = BigInt::from(156);
    generators[6][0] = BigInt::from(-28826);
    generators[6][1] = BigInt::from(-148);
    generators[6][2] = BigInt::from(220);
    generators[6][3] = BigInt::from(-122);
    generators[7][0] = BigInt::from(586);
    generators[7][1] = BigInt::from(-28426);
    generators[7][2] = BigInt::from(278);
    generators[7][3] = BigInt::from(218);

    let mut expected_cmp = MatrixUtils::zero::<BigInt>();
    expected_cmp[0][0] = BigInt::from(2321156);
    expected_cmp[1][1] = BigInt::from(2321156);
    expected_cmp[0][2] = BigInt::from(620252);
    expected_cmp[1][2] = BigInt::from(365058);
    expected_cmp[2][2] = BigInt::from(2);
    expected_cmp[0][3] = BigInt::from(1956098);
    expected_cmp[1][3] = BigInt::from(620252);
    expected_cmp[3][3] = BigInt::from(2);

    let det = BigInt::parse_bytes(b"21551060705344", 10).unwrap();
    let mut hnf = MatrixUtils::zero::<BigInt>();
    
    mat_4xn_hnf_mod_core(&mut hnf, &generators, &det);
    
    assert!(MatrixUtils::equal(&expected_cmp, &hnf), "GitHub #38 exact vector test failed");
}

#[test]
fn test_mat_4xn_hnf_mod_core_fuzzing() {
    let mut rng = StdRng::seed_from_u64(123456789);

    for iteration in 0..250 {
        let n_gens = rng.random_range(4..12);

        // 1. Generate a mathematically valid SQISign HNF matrix
        let mut expected_hnf = MatrixUtils::zero::<BigInt>();
        
        for i in 0..4 {
            expected_hnf[i][i] = BigInt::from_i32(rng.random_range(1..25));
        }
        
        for i in 0..4 {
            for j in (i + 1)..4 {
                let d_val = i32::try_from(&expected_hnf[i][i]).unwrap_or(1);
                expected_hnf[i][j] = BigInt::from_i32(rng.random_range(0..d_val));
            }
        }

        assert_is_hnf(&expected_hnf);

        let mut det = b_one();
        for i in 0..4 {
            det = det * expected_hnf[i][i].clone();
        }

        let mut gens = Vec::with_capacity(n_gens);
        for _ in 0..n_gens {
            gens.push([b_zero(), b_zero(), b_zero(), b_zero()]);
        }
        for i in 0..4 {
            for j in 0..4 {
                gens[i][j] = expected_hnf[j][i].clone();
            }
        }

        let n_ops = 300;
        for _ in 0..n_ops {
            let op = rng.random_range(0..3);
            if op == 0 {
                let r1 = rng.random_range(0..n_gens);
                let r2 = rng.random_range(0..n_gens);
                if r1 != r2 {
                    let scalar = BigInt::from_i32(rng.random_range(-5..5));
                    for c in 0..4 {
                        gens[r2][c] = gens[r2][c].clone() + scalar.clone() * gens[r1][c].clone();
                    }
                }
            } else if op == 1 {
                let r1 = rng.random_range(0..n_gens);
                let r2 = rng.random_range(0..n_gens);
                if r1 != r2 {
                    gens.swap(r1, r2);
                }
            } else {
                let r1 = rng.random_range(0..n_gens);
                for c in 0..4 {
                    gens[r1][c] = -gens[r1][c].clone();
                }
            }
        }

        for i in 0..n_gens {
            for c in 0..4 {
                let scalar = BigInt::from_i32(rng.random_range(-3..4));
                gens[i][c] = gens[i][c].clone() + scalar * det.clone();
            }
        }

        let mut computed_hnf = MatrixUtils::zero::<BigInt>();
        mat_4xn_hnf_mod_core(&mut computed_hnf, &gens, &det);

        if !MatrixUtils::equal(&expected_hnf, &computed_hnf) {
            println!("\n[!] FUZZING FAILED ON ITERATION {}", iteration);
            println!("Determinant (Modulo): {}", det);
            println!("Expected HNF Matrix:");
            for i in 0..4 { println!("  {:?}", expected_hnf[i]); }
            println!("Computed HNF Matrix:");
            for i in 0..4 { println!("  {:?}", computed_hnf[i]); }
            panic!("mat_4xn_hnf_mod_core failed to reconstruct the original lattice.");
        }
    }
}
