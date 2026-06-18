// src/quaternion/lattice.rs

use crate::quaternion::algebra::{BigIntAlg, QuatAlg, QuatElem};
use crate::quaternion::hnf::{mat_4x4_inv_with_det_as_denom, mat_4xn_hnf_mod_core};
use crate::quaternion::matrix::{Mat4x4, MatrixUtils, Vec4};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatLattice<T: BigIntAlg> {
    /// TODO: use a Vec4 of quaternions instead. maybe?
    /// would be cool to have "ZZ quaternion" and "QQ quaternion = ZZ quaternion + denom"
    /// so that quat lattice is exactly "4 × ZZ quaternion + denom"
    pub basis: Mat4x4<T>,
    pub denom: T,
}

impl<T: BigIntAlg> QuatLattice<T> {
    pub fn new(basis: Mat4x4<T>, denom: T) -> Self {
        debug_assert!(!denom.is_zero(), "Lattice denominator cannot be zero");
        Self { basis, denom }
    }

    pub fn zero() -> Self {
        Self {
            basis: MatrixUtils::zero(),
            denom: T::one(),
        }
    }

    #[inline]
    fn coord_mul(a: &[T; 4], b: &[T; 4], p: &T) -> [T; 4] {
        let c0 = a[0].clone() * b[0].clone() - a[1].clone() * b[1].clone()
            - p.clone() * (a[2].clone() * b[2].clone() + a[3].clone() * b[3].clone());

        let c1 = p.clone() * (a[2].clone() * b[3].clone() - a[3].clone() * b[2].clone())
            + a[0].clone() * b[1].clone() + a[1].clone() * b[0].clone();

        let c2 = a[0].clone() * b[2].clone() + a[2].clone() * b[0].clone()
            - a[1].clone() * b[3].clone() + a[3].clone() * b[1].clone();

        let c3 = a[0].clone() * b[3].clone() + a[3].clone() * b[0].clone()
            - a[2].clone() * b[1].clone() + a[1].clone() * b[2].clone();

        [c0, c1, c2, c3]
    }

    pub fn reduce_denom(&mut self) {
        let mut g = self.denom.clone();
        for i in 0..4 {
            for j in 0..4 {
                g = g.gcd(&self.basis[i][j]);
            }
        }

        let mut sign = T::one();
        if self.denom < T::zero() {
            sign = T::from_i32(-1);
        }

        let divisor = g * sign;
        for i in 0..4 {
            for j in 0..4 {
                self.basis[i][j] = self.basis[i][j].clone() / divisor.clone();
            }
        }
        self.denom = self.denom.clone() / divisor;
    }

    pub fn hnf(&mut self) {
        // following quat_lattice_hnf of sqisign… maybe it would be more
        // convenient to have either hnf or lattice transposed to avoid
        // having to do this?
        let mod_val = mat_4x4_inv_with_det_as_denom(None, &self.basis).abs();

        let mut generators = Vec::with_capacity(4);
        for j in 0..4 {
            let mut row = [T::zero(), T::zero(), T::zero(), T::zero()];
            for i in 0..4 {
                row[i] = self.basis[i][j].clone();
            }
            generators.push(row);
        }

        mat_4xn_hnf_mod_core(&mut self.basis, &generators, &mod_val);
        self.reduce_denom();
    }

    pub fn equal(lat1: &Self, lat2: &Self) -> bool {
        let mut a = lat1.clone();
        let mut b = lat2.clone();
        a.reduce_denom();
        b.reduce_denom();
        a.denom = a.denom.abs();
        b.denom = b.denom.abs();
        a.hnf();
        b.hnf();
        // TODO: That's two reduce_denom, probably useless

        (a.denom == b.denom) && MatrixUtils::equal(&a.basis, &b.basis)
    }

    pub fn inclusion(sublat: &Self, overlat: &Self) -> bool {
        let sum = Self::add(overlat, sublat);
        Self::equal(&sum, overlat)
    }

    pub fn conjugate_without_hnf(&self) -> Self {
        let mut conj = self.clone();
        for row in 1..4 {
            for col in 0..4 {
                conj.basis[row][col] = -conj.basis[row][col].clone();
            }
        }
        conj
    }

    pub fn dual_without_hnf(&self) -> Self {
        let mut inv = MatrixUtils::zero::<T>();
        let det = mat_4x4_inv_with_det_as_denom(Some(&mut inv), &self.basis);

        let mut inv_t = MatrixUtils::zero::<T>();
        for i in 0..4 {
            for j in 0..4 {
                inv_t[i][j] = inv[j][i].clone();
            }
        }

        let mut basis = MatrixUtils::zero::<T>();
        for i in 0..4 {
            for j in 0..4 {
                basis[i][j] = self.denom.clone() * inv_t[i][j].clone();
            }
        }

        Self { basis, denom: det }
    }

    pub fn add(lat1: &Self, lat2: &Self) -> Self {
        let mut generators = Vec::with_capacity(8);

        // tmp1 = lat1.denom * lat2.basis
        let mut tmp1 = MatrixUtils::zero::<T>();
        for i in 0..4 {
            for j in 0..4 {
                tmp1[i][j] = lat1.denom.clone() * lat2.basis[i][j].clone();
            }
        }
        for j in 0..4 {
            let mut row = [T::zero(), T::zero(), T::zero(), T::zero()];
            for i in 0..4 {
                row[i] = tmp1[i][j].clone();
            }
            generators.push(row);
        }
        let det1 = mat_4x4_inv_with_det_as_denom(None, &tmp1);

        // tmp2 = lat2.denom * lat1.basis
        let mut tmp2 = MatrixUtils::zero::<T>();
        for i in 0..4 {
            for j in 0..4 {
                tmp2[i][j] = lat2.denom.clone() * lat1.basis[i][j].clone();
            }
        }
        for j in 0..4 {
            let mut row = [T::zero(), T::zero(), T::zero(), T::zero()];
            for i in 0..4 {
                row[i] = tmp2[i][j].clone();
            }
            generators.push(row);
        }
        let det2 = mat_4x4_inv_with_det_as_denom(None, &tmp2);

        debug_assert!(!det1.is_zero(), "Lattice 1 det is zero in add()");
        debug_assert!(!det2.is_zero(), "Lattice 2 det is zero in add()");

        let detprod = det1.gcd(&det2);
        let mut res = QuatLattice::zero();

        mat_4xn_hnf_mod_core(&mut res.basis, &generators, &detprod);
        res.denom = lat1.denom.clone() * lat2.denom.clone();
        res.reduce_denom();

        res
    }

    pub fn intersect(lat1: &Self, lat2: &Self) -> Self {
        let dual1 = lat1.dual_without_hnf();
        let dual2 = lat2.dual_without_hnf();
        let dual_res = Self::add(&dual1, &dual2);

        let mut res = dual_res.dual_without_hnf();
        res.hnf(); // Maybe not necessary, let callers do that if they want TODO choose design
        res
    }

    fn mat_alg_coord_mul_without_hnf(
        lat: &Mat4x4<T>,
        coord: &Vec4<T>,
        alg: &QuatAlg<T>,
    ) -> Mat4x4<T> {
        let mut prod = MatrixUtils::zero::<T>();
        for i in 0..4 {
            let mut a = [T::zero(), T::zero(), T::zero(), T::zero()];
            for j in 0..4 {
                a[j] = lat[j][i].clone();
            }
            // TODO: this is stupid
            let p = Self::coord_mul(&a, coord, &alg.p);
            for j in 0..4 {
                prod[j][i] = p[j].clone();
            }
        }
        prod
    }

    pub fn alg_elem_mul(lat: &Self, elem: &QuatElem<T>, alg: &QuatAlg<T>) -> Self {
        // TODO: this is stupid
        let mut res = QuatLattice::zero();
        res.basis = Self::mat_alg_coord_mul_without_hnf(&lat.basis, &elem.coord, alg);
        res.denom = lat.denom.clone() * elem.denom.clone();
        res.hnf();
        res
    }

    pub fn mul(lat1: &Self, lat2: &Self, alg: &QuatAlg<T>) -> Self {
        let mut generators = Vec::with_capacity(16);
        let mut detmat = MatrixUtils::zero::<T>();

        for k in 0..4 {
            let mut elem1 = [T::zero(), T::zero(), T::zero(), T::zero()];
            for r in 0..4 {
                elem1[r] = lat1.basis[r][k].clone();
            }

            for i in 0..4 {
                let mut elem2 = [T::zero(), T::zero(), T::zero(), T::zero()];
                for r in 0..4 {
                    elem2[r] = lat2.basis[r][i].clone();
                }

                let elem_res = Self::coord_mul(&elem1, &elem2, &alg.p);

                for j in 0..4 {
                    if k == 0 {
                        detmat[i][j] = elem_res[j].clone();
                    }
                }
                generators.push(elem_res);
            }
        }

        let det = mat_4x4_inv_with_det_as_denom(None, &detmat).abs();
        let mut res = QuatLattice::zero();
        mat_4xn_hnf_mod_core(&mut res.basis, &generators, &det);
        res.denom = lat1.denom.clone() * lat2.denom.clone();
        res.reduce_denom();

        res
    }

    pub fn contains(&self, x: &QuatElem<T>) -> Option<Vec4<T>> {
        let mut inv = MatrixUtils::zero::<T>();
        let det = mat_4x4_inv_with_det_as_denom(Some(&mut inv), &self.basis);
        debug_assert!(!det.is_zero(), "Lattice basis determinant is zero in contains()");

        let mut work_coord = [T::zero(), T::zero(), T::zero(), T::zero()];
        for i in 0..4 {
            let mut sum = T::zero();
            for j in 0..4 {
                sum = sum + inv[i][j].clone() * x.coord[j].clone();
            }
            work_coord[i] = sum;
        }

        for i in 0..4 {
            work_coord[i] = work_coord[i].clone() * self.denom.clone();
        }

        let prod = x.denom.clone() * det;
        let mut divisible = true;

        for i in 0..4 {
            if !(work_coord[i].clone() % prod.clone()).is_zero() {
                divisible = false;
            }
            work_coord[i] = work_coord[i].clone() / prod.clone();
        }

        if divisible {
            Some(work_coord)
        } else {
            None
        }
    }

    pub fn index(sublat: &Self, overlat: &Self) -> T {
        let det_sub = mat_4x4_inv_with_det_as_denom(None, &sublat.basis);
        let mut tmp_over = overlat.denom.clone() * overlat.denom.clone();
        tmp_over = tmp_over.clone() * tmp_over.clone();
        let num = det_sub * tmp_over;

        let mut tmp_sub = sublat.denom.clone() * sublat.denom.clone();
        tmp_sub = tmp_sub.clone() * tmp_sub.clone();
        let det_over = mat_4x4_inv_with_det_as_denom(None, &overlat.basis);
        let den = tmp_sub * det_over;

        let index = num / den;
        index.abs()
    }

    pub fn gram(&self, alg: &QuatAlg<T>) -> Mat4x4<T> {
        let mut g_mat = MatrixUtils::zero::<T>();
        let two = T::from_i32(2);

        for i in 0..4 {
            for j in 0..=i {
                let mut sum = T::zero();
                for k in 0..4 {
                    let mut tmp = self.basis[k][i].clone() * self.basis[k][j].clone();
                    if k >= 2 {
                        tmp = tmp * alg.p.clone();
                    }
                    sum = sum + tmp;
                }
                g_mat[i][j] = sum * two.clone();
            }
        }

        for i in 0..4 {
            for j in (i + 1)..4 {
                g_mat[i][j] = g_mat[j][i].clone();
            }
        }
        g_mat
    }
}
