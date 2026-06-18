// src/quaternion/hnf.rs

use crate::quaternion::algebra::BigIntAlg;
use crate::quaternion::matrix::{Mat4x4, MatrixUtils, Vec4};

/// Euclidean modulo that always returns a positive remainder in [0, modulo - 1].
/// TODO: is it really necessary? maybe force this to be true in BigIntAlg
/// maybe better architecture to have these directly in BigIntAlg traits as well
fn euclidean_mod<T: BigIntAlg>(x: &T, modulo: &T) -> T {
    let mut r = x.clone() % modulo.clone();
    if r < T::zero() {
        r = r + modulo.clone();
    }
    r
}

/// Modulo in [1, modulo] instead
/// maybe better architecture to have these directly in BigIntAlg traits as well
pub fn mod_not_zero<T: BigIntAlg>(x: &T, modulo: &T) -> T {
    let mut r = euclidean_mod(x, modulo);
    if r.is_zero() {
        r = modulo.clone();
    }
    r
}

/// Centered modulo (results in range (-mod/2, mod/2])
pub fn centered_mod<T: BigIntAlg>(a: &T, modulo: &T) -> T {
    let two = T::from_i32(2);
    let d = modulo.clone() / two; // Floor division for integers

    let tmp = mod_not_zero(a, modulo);

    if tmp > d {
        tmp - modulo.clone()
    } else {
        tmp
    }
}

/// Extended GCD except the first coefficient (u) is not zero in the result
pub fn xgcd_with_u_not_0<T: BigIntAlg>(a: &T, b: &T) -> (T, T, T) {
    let (d, mut u, mut v) = a.xgcd(b);

    if u.is_zero() {
        let t1 = b.clone() / d.clone();
        u = u + t1;

        let t2 = a.clone() / d.clone();
        v = v - t2;
    }

    (d, u, v)
}

/// 3x3 matrix determinant TODO: maybe move that to matrix.rs
fn det_3x3<T: BigIntAlg>(m: &[[T; 3]; 3]) -> T {
    let term1 = m[0][0].clone() * (m[1][1].clone() * m[2][2].clone() - m[1][2].clone() * m[2][1].clone());
    let term2 = m[0][1].clone() * (m[1][0].clone() * m[2][2].clone() - m[1][2].clone() * m[2][0].clone());
    let term3 = m[0][2].clone() * (m[1][0].clone() * m[2][1].clone() - m[1][1].clone() * m[2][0].clone());
    term1 - term2 + term3
}

/// Computes the adjoint matrix (if requested) and the determinant.
/// TODO: maybe move that to matrix.rs
/// TODO: this is naive way - sqisign uses better algorithm
pub fn mat_4x4_inv_with_det_as_denom<T: BigIntAlg>(
    inv: Option<&mut Mat4x4<T>>,
    mat: &Mat4x4<T>,
) -> T {
    let mut adj = MatrixUtils::zero::<T>();

    // Calculate cofactors
    for i in 0..4 {
        for j in 0..4 {
            let mut minor = [
                [T::zero(), T::zero(), T::zero()],
                [T::zero(), T::zero(), T::zero()],
                [T::zero(), T::zero(), T::zero()],
            ];

            let mut mi = 0;
            for row in 0..4 {
                if row == i { continue; }
                let mut mj = 0;
                for col in 0..4 {
                    if col == j { continue; }
                    minor[mi][mj] = mat[row][col].clone();
                    mj += 1;
                }
                mi += 1;
            }

            let mut c = det_3x3(&minor);
            if (i + j) % 2 != 0 {
                c = -c;
            }

            adj[j][i] = c;
        }
    }

    let mut det = T::zero();
    for k in 0..4 {
        det = det + mat[0][k].clone() * adj[k][0].clone();
    }

    if let Some(inv_mat) = inv {
        *inv_mat = adj;
    }

    det
}

/// Hermite Normal Form Modulo operation.
pub fn mat_4xn_hnf_mod_core<T: BigIntAlg>(
    hnf: &mut Mat4x4<T>,
    generators: &[Vec4<T>],
    modulo: &T,
) {
    let n = generators.len();
    debug_assert!(n >= 4, "HNF modulo requires at least 4 generators");

    let mut a = Vec::with_capacity(n);
    for g in generators {
        a.push(g.clone());
    }

    let mut w = MatrixUtils::zero::<T>();
    let mut m = modulo.clone();
    if m < T::zero() {
        m = -m;
    }

    let mut i: isize = 3;
    let mut k: isize = (n as isize) - 1;
    let mut j: isize = (n as isize) - 1;

    while i >= 0 {
        let ui = i as usize;
        let uk = k as usize;

        while j > 0 {
            j -= 1;
            let uj = j as usize;

            if !a[uj][ui].is_zero() {
                let (d, u, v) = xgcd_with_u_not_0(&a[uk][ui], &a[uj][ui]);

                // c = u*a[k] + v*a[j]
                let mut c = [T::zero(), T::zero(), T::zero(), T::zero()];
                for x in 0..4 {
                    c[x] = u.clone() * a[uk][x].clone() + v.clone() * a[uj][x].clone();
                }

                let coeff_1 = a[uk][ui].clone() / d.clone();
                let coeff_2 = -(a[uj][ui].clone() / d.clone());

                // a[j] = coeff_1*a[j] + coeff_2*a[k] mod m
                for x in 0..4 {
                    let sum = coeff_1.clone() * a[uj][x].clone() + coeff_2.clone() * a[uk][x].clone();
                    a[uj][x] = centered_mod(&sum, &m);
                }

                // a[k] = c mod m TODO: check is it centered mod? I think
                for x in 0..4 {
                    a[uk][x] = centered_mod(&c[x], &m);
                }
            }
        }

        let (d, u, _v) = xgcd_with_u_not_0(&a[uk][ui], &m);

        // w[i] = u * a[k] mod m (positive remainder mod)
        for x in 0..4 {
            let mut prod = u.clone() * a[uk][x].clone();
            // TODO: use mod not zero
            prod = prod % m.clone();
            if prod < T::zero() {
                prod = prod + m.clone();
            }
            w[ui][x] = prod;
        }

        // TODO: quite sure this is useless
        if w[ui][ui].is_zero() {
            w[ui][ui] = m.clone();
        }

        for h in (ui + 1)..4 {
            // q = floor(w[h][i] / w[i][i])
            let mut q = w[h][ui].clone() / w[ui][ui].clone();
            let rem = w[h][ui].clone() % w[ui][ui].clone();
            if rem < T::zero() {
                q = q - T::from_i32(1); // floor adjustment
            }
            q = -q;

            for x in 0..4 {
                w[h][x] = w[h][x].clone() + q.clone() * w[ui][x].clone();
            }
        }

        m = m / d;

        if i != 0 {
            k -= 1;
            i -= 1;
            j = k;
            if a[k as usize][i as usize].is_zero() {
                a[k as usize][i as usize] = m.clone();
            }
        } else {
            k -= 1;
            i -= 1;
            j = k;
        }
    }

    // transpose? I guess
    for x in 0..4 {
        for y in 0..4 {
            hnf[x][y] = w[y][x].clone();
        }
    }
}
