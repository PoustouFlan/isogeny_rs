use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
use core::fmt::Debug;

/// Generic interface to allow switching between num_bigint, rug, or crypto_bigint.
pub trait BigIntAlg:
    Clone
    + Debug
    + PartialEq
    + Eq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + Neg<Output = Self>
    + Sized
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_i32(val: i32) -> Self;

    fn gcd(&self, other: &Self) -> Self;

    fn abs(&self) -> Self;
    fn is_zero(&self) -> bool;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatAlg<T: BigIntAlg> {
    pub p: T,
}

impl<T: BigIntAlg> QuatAlg<T> {
    pub fn new(p: T) -> Self {
        Self { p }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatElem<T: BigIntAlg> {
    pub denom: T,
    pub coord: [T; 4],
}

impl<T: BigIntAlg> QuatElem<T> {
    pub fn new_i32(denom: i32, c0: i32, c1: i32, c2: i32, c3: i32) -> Self {
        debug_assert!(denom != 0, "QuatElem initialized with zero denominator");
        Self {
            denom: T::from_i32(denom),
            coord: [
                T::from_i32(c0),
                T::from_i32(c1),
                T::from_i32(c2),
                T::from_i32(c3),
            ],
        }
    }

    pub fn from_scalar(numerator: T, denominator: T) -> Self {
        debug_assert!(!denominator.is_zero(), "QuatElem initialized with zero denominator");
        Self {
            denom: denominator,
            coord: [numerator, T::zero(), T::zero(), T::zero()],
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coord[0].is_zero()
            && self.coord[1].is_zero()
            && self.coord[2].is_zero()
            && self.coord[3].is_zero()
    }

    pub fn conj(&self) -> Self {
        Self {
            denom: self.denom.clone(),
            coord: [
                self.coord[0].clone(),
                -self.coord[1].clone(),
                -self.coord[2].clone(),
                -self.coord[3].clone(),
            ],
        }
    }

    /// Aligns denominators of two elements. Returns (scaled_a, scaled_b).
    fn equal_denom(a: &Self, b: &Self) -> (Self, Self) {
        debug_assert!(!a.denom.is_zero(), "a.denom is zero in equal_denom");
        debug_assert!(!b.denom.is_zero(), "b.denom is zero in equal_denom");

        let gcd = a.denom.gcd(&b.denom);
        let a_mult = b.denom.clone() / gcd.clone();
        let b_mult = a.denom.clone() / gcd;

        let mut res_a = a.clone();
        let mut res_b = b.clone();

        for i in 0..4 {
            res_a.coord[i] = res_a.coord[i].clone() * a_mult.clone();
            res_b.coord[i] = res_b.coord[i].clone() * b_mult.clone();
        }

        let common_denom = a.denom.clone() * a_mult;
        res_a.denom = common_denom.clone();
        res_b.denom = common_denom;

        (res_a, res_b)
    }

    pub fn add(&self, other: &Self) -> Self {
        let (mut a_scaled, b_scaled) = Self::equal_denom(self, other);
        for i in 0..4 {
            a_scaled.coord[i] = a_scaled.coord[i].clone() + b_scaled.coord[i].clone();
        }
        a_scaled
    }

    pub fn sub(&self, other: &Self) -> Self {
        let (mut a_scaled, b_scaled) = Self::equal_denom(self, other);
        for i in 0..4 {
            a_scaled.coord[i] = a_scaled.coord[i].clone() - b_scaled.coord[i].clone();
        }
        a_scaled
    }

    #[inline]
    fn coord_mul(a: &[T; 4], b: &[T; 4], p: &T) -> [T; 4] {
        // c0 = a0b0 - a1b1 - p(a2b2 + a3b3)
        let mut c0 = a[0].clone() * b[0].clone() - a[1].clone() * b[1].clone();
        c0 = c0 - p.clone() * (a[2].clone() * b[2].clone() + a[3].clone() * b[3].clone());

        // c1 = p(a2b3 - a3b2) + a0b1 + a1b0
        let mut c1 = p.clone() * (a[2].clone() * b[3].clone() - a[3].clone() * b[2].clone());
        c1 = c1 + a[0].clone() * b[1].clone() + a[1].clone() * b[0].clone();

        // c2 = a0b2 + a2b0 - a1b3 + a3b1
        let c2 = a[0].clone() * b[2].clone() + a[2].clone() * b[0].clone()
            - a[1].clone() * b[3].clone() + a[3].clone() * b[1].clone();

        // c3 = a0b3 + a3b0 - a2b1 + a1b2
        let c3 = a[0].clone() * b[3].clone() + a[3].clone() * b[0].clone()
            - a[2].clone() * b[1].clone() + a[1].clone() * b[2].clone();

        [c0, c1, c2, c3]
    }

    pub fn mul(&self, other: &Self, alg: &QuatAlg<T>) -> Self {
        let denom = self.denom.clone() * other.denom.clone();
        let coord = Self::coord_mul(&self.coord, &other.coord, &alg.p);

        Self { denom, coord }
    }

    /// Returns (numerator, denominator) of the reduced norm.
    pub fn norm(&self, alg: &QuatAlg<T>) -> (T, T) {
        let conj_x = self.conj();
        let norm_elem = self.mul(&conj_x, alg);

        let num = norm_elem.coord[0].clone();
        let denom = norm_elem.denom.clone();

        let g = num.gcd(&denom);
        let res_num = (num / g.clone()).abs();
        let res_denom = (denom / g).abs();

        (res_num, res_denom)
    }

    pub fn normalize(&mut self) {
        debug_assert!(!self.denom.is_zero(), "Denominator is zero during normalize");

        let mut g = self.coord[0].clone();
        for i in 1..4 {
            g = g.gcd(&self.coord[i]);
        }
        g = g.gcd(&self.denom);

        let mut sign = T::one();

        // TODO: make that constant time
        if self.denom < T::zero() {
            sign = T::from_i32(-1);
        }

        let divisor = g * sign;

        for i in 0..4 {
            self.coord[i] = self.coord[i].clone() / divisor.clone();
        }
        self.denom = self.denom.clone() / divisor;
    }
}
