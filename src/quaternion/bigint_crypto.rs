use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
use crypto_bigint::{Int, NonZero};
use crate::quaternion::algebra::BigIntAlg;

/// Wrapper around crypto_bigint::Int to satisfy the standard Rust 
/// operators required by the generic `BigIntAlg` trait.
/// Newtype pattern should be zero-cost abstraction.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CryptoInt<const LIMBS: usize>(pub Int<LIMBS>);

impl<const LIMBS: usize> Add for CryptoInt<LIMBS> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        CryptoInt(self.0 + rhs.0)
    }
}

impl<const LIMBS: usize> Sub for CryptoInt<LIMBS> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        CryptoInt(self.0 - rhs.0)
    }
}

impl<const LIMBS: usize> Mul for CryptoInt<LIMBS> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        CryptoInt(self.0 * rhs.0)
    }
}

impl<const LIMBS: usize> Div for CryptoInt<LIMBS> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let nz: Option<NonZero<Int<LIMBS>>> = NonZero::new(rhs.0).into();
        let res: Option<Int<LIMBS>> = (self.0 / nz.expect("Division by zero in CryptoInt::div")).into();
        CryptoInt(res.expect("Division overflowed"))
    }
}

impl<const LIMBS: usize> Rem for CryptoInt<LIMBS> {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self::Output {
        let nz: Option<NonZero<Int<LIMBS>>> = NonZero::new(rhs.0).into();
        let res: Option<Int<LIMBS>> = (self.0 % nz.expect("Division by zero in CryptoInt::rem")).into();
        CryptoInt(res.expect("Remainder overflowed"))
    }
}

impl<const LIMBS: usize> Neg for CryptoInt<LIMBS> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        // Handled via zero subtraction to emulate 2's complement negation
        CryptoInt(Int::ZERO - self.0)
    }
}

impl<const LIMBS: usize> BigIntAlg for CryptoInt<LIMBS> {
    fn zero() -> Self {
        CryptoInt(Int::ZERO)
    }

    fn one() -> Self {
        CryptoInt(Int::ONE)
    }

    fn from_i32(val: i32) -> Self {
        CryptoInt(Int::from_i32(val))
    }

    fn gcd(&self, other: &Self) -> Self {
        let a_abs = self.0.abs();
        let b_abs = other.0.abs();

        // Note: unsigned GCD
        let g_uint = a_abs.gcd_vartime(&b_abs);
        // Convert Uint back to Int by passing the limbs array
        CryptoInt(Int::new(g_uint.into()))
    }

    fn abs(&self) -> Self {
        let uint_abs = self.0.abs();
        // Convert Uint back to Int by passing the limbs array
        CryptoInt(Int::new(uint_abs.into()))
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero().into()
    }
}
