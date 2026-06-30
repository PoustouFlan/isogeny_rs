use core::ops::{Add, Mul};
use crate::quaternion::algebra::{BigIntAlg, QuatConfig, IntQuat, RatQuat};
use crate::quaternion::hnf::quat_hnf_mod_core;
use crate::quaternion::matrix::MatrixUtils;

#[derive(Debug)]
pub struct QuatLattice<T: BigIntAlg, P: QuatConfig<T>> {
    pub generators: [IntQuat<T, P>; 4],
    pub denom: T,
}

impl<T: BigIntAlg, P: QuatConfig<T>> Clone for QuatLattice<T, P> {
    fn clone(&self) -> Self {
        Self {
            generators: self.generators.clone(),
            denom: self.denom.clone(),
        }
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> PartialEq for QuatLattice<T, P> {
    fn eq(&self, other: &Self) -> bool {
        self.denom == other.denom && self.generators == other.generators
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> Eq for QuatLattice<T, P> {}

impl<T: BigIntAlg, P: QuatConfig<T>> QuatLattice<T, P> {
    pub fn zero() -> Self {
        Self {
            generators: [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()],
            denom: T::one(),
        }
    }

    pub fn reduce_denom(&mut self) {
        let mut g = self.denom.clone();
        for i in 0..4 {
            g = g.gcd(&self.generators[i].coords_gcd());
        }

        let mut sign = T::one();
        if self.denom < T::zero() {
            sign = T::from_i32(-1);
        }

        let divisor = g * sign;
        for i in 0..4 {
            self.generators[i] = &self.generators[i] / &divisor;
        }
        self.denom = self.denom.clone() / divisor;
    }

    pub fn hnf(&mut self) {
        let mod_val = MatrixUtils::mat_4x4_inv_with_det_as_denom(None, &self.generators).abs();
        self.generators = quat_hnf_mod_core(&mut self.generators, &mod_val);
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

        a.denom == b.denom && a.generators == b.generators
    }

    pub fn inclusion(sublat: &Self, overlat: &Self) -> bool {
        let sum = Self::add_lazy(overlat, sublat);
        Self::equal(&sum, overlat)
    }

    pub fn conjugate_without_hnf(&self) -> Self {
        Self {
            generators: [
                self.generators[0].conj(),
                self.generators[1].conj(),
                self.generators[2].conj(),
                self.generators[3].conj(),
            ],
            denom: self.denom.clone(),
        }
    }

    pub fn dual_without_hnf(&self) -> Self {
        let mut inv = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];
        let det = MatrixUtils::mat_4x4_inv_with_det_as_denom(Some(&mut inv), &self.generators);

        for i in 0..4 {
            inv[i] = &inv[i] * &self.denom;
        }

        Self { generators: inv, denom: det }
    }

    pub fn add_lazy(lat1: &Self, lat2: &Self) -> Self {
        let mut generators = Vec::with_capacity(8);
        for i in 0..4 { generators.push(&lat2.generators[i] * &lat1.denom); }
        for i in 0..4 { generators.push(&lat1.generators[i] * &lat2.denom); }

        // Since we dropped Mat4x4, we rebuild the det matrices here for the mod_val
        let mut mat1 = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];
        let mut mat2 = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];
        for i in 0..4 { mat1[i] = generators[i].clone(); mat2[i] = generators[i+4].clone(); }

        let det1 = MatrixUtils::mat_4x4_inv_with_det_as_denom(None, &mat1);
    let det2 = MatrixUtils::mat_4x4_inv_with_det_as_denom(None, &mat2);

        let mut res = Self {
            generators: quat_hnf_mod_core(&mut generators, &det1.gcd(&det2)),
            denom: lat1.denom.clone() * lat2.denom.clone(),
        };
        res.reduce_denom();
        res
    }

    pub fn intersect(lat1: &Self, lat2: &Self) -> Self {
        let dual1 = lat1.dual_without_hnf();
        let dual2 = lat2.dual_without_hnf();
        let dual_res = Self::add_lazy(&dual1, &dual2);

        let mut res = dual_res.dual_without_hnf();
        res.hnf(); 
        res
    }

    pub fn alg_elem_mul(lat: &Self, elem: &RatQuat<T, P>) -> Self {
        let mut generators = Vec::with_capacity(4);
        for i in 0..4 { generators.push(&lat.generators[i] * &elem.num); }

        let mut mat = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];
        for i in 0..4 { mat[i] = generators[i].clone(); }

        let mod_val = MatrixUtils::mat_4x4_inv_with_det_as_denom(None, &mat).abs();
        let mut res = Self {
            generators: quat_hnf_mod_core(&mut generators, &mod_val),
            denom: lat.denom.clone() * elem.denom.clone(),
        };
        res.reduce_denom();
        res
    }

    pub fn mul_lazy(lat1: &Self, lat2: &Self) -> Self {
        let mut generators = Vec::with_capacity(16);
        let mut detmat = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];

        for k in 0..4 {
            for i in 0..4 {
                let elem_res = &lat1.generators[k] * &lat2.generators[i];
                if k == 0 { detmat[i] = elem_res.clone(); }
                generators.push(elem_res);
            }
        }

        let mod_val = MatrixUtils::mat_4x4_inv_with_det_as_denom(None, &detmat).abs();
        let mut res = Self {
            generators: quat_hnf_mod_core(&mut generators, &mod_val),
            denom: lat1.denom.clone() * lat2.denom.clone(),
        };
        res.reduce_denom();
        res
    }
}

// Eager operations bound to standard Rust operators
impl<T: BigIntAlg, P: QuatConfig<T>> Add<QuatLattice<T, P>> for QuatLattice<T, P> {
    type Output = QuatLattice<T, P>;
    fn add(self, rhs: QuatLattice<T, P>) -> Self::Output {
        let mut res = QuatLattice::add_lazy(&self, &rhs);
        res.hnf();
        res
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Add<&'b QuatLattice<T, P>> for &'a QuatLattice<T, P> {
    type Output = QuatLattice<T, P>;
    fn add(self, rhs: &'b QuatLattice<T, P>) -> Self::Output {
        let mut res = QuatLattice::add_lazy(self, rhs);
        res.hnf();
        res
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> Mul<QuatLattice<T, P>> for QuatLattice<T, P> {
    type Output = QuatLattice<T, P>;
    fn mul(self, rhs: QuatLattice<T, P>) -> Self::Output {
        let mut res = QuatLattice::mul_lazy(&self, &rhs);
        res.hnf();
        res
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Mul<&'b QuatLattice<T, P>> for &'a QuatLattice<T, P> {
    type Output = QuatLattice<T, P>;
    fn mul(self, rhs: &'b QuatLattice<T, P>) -> Self::Output {
        let mut res = QuatLattice::mul_lazy(self, rhs);
        res.hnf();
        res
    }
}
