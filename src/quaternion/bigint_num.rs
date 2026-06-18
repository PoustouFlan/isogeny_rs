use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero, One};
use crate::quaternion::algebra::BigIntAlg;

/// Implementation of BigIntAlg that uses num_bigint backend
impl BigIntAlg for BigInt {
    fn zero() -> Self {
        <BigInt as Zero>::zero()
    }

    fn one() -> Self {
        <BigInt as One>::one()
    }

    fn from_i32(val: i32) -> Self {
        BigInt::from(val)
    }

    fn gcd(&self, other: &Self) -> Self {
        self.extended_gcd(other).gcd
    }

    fn xgcd(&self, other: &Self) -> (Self, Self, Self) {
        let res = self.extended_gcd(other);
        (res.gcd, res.x, res.y)
    }

    fn abs(&self) -> Self {
        <BigInt as Signed>::abs(self)
    }

    fn is_zero(&self) -> bool {
        <BigInt as Zero>::is_zero(self)
    }
}
