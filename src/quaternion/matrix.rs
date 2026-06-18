use crate::quaternion::algebra::BigIntAlg;

pub type Vec4<T> = [T; 4];
pub type Mat4x4<T> = [Vec4<T>; 4];

/// Helper functions for 4x4 matrices over BigIntAlg
pub struct MatrixUtils;

impl MatrixUtils {
    /// Creates a 4x4 matrix initialized to zero
    pub fn zero<T: BigIntAlg>() -> Mat4x4<T> {
        [
            [T::zero(), T::zero(), T::zero(), T::zero()],
            [T::zero(), T::zero(), T::zero(), T::zero()],
            [T::zero(), T::zero(), T::zero(), T::zero()],
            [T::zero(), T::zero(), T::zero(), T::zero()],
        ]
    }

    /// Creates a 4x4 identity matrix
    pub fn identity<T: BigIntAlg>() -> Mat4x4<T> {
        let mut mat = Self::zero();
        mat[0][0] = T::one();
        mat[1][1] = T::one();
        mat[2][2] = T::one();
        mat[3][3] = T::one();
        mat
    }

    /// Checks if two matrices are exactly equal
    pub fn equal<T: BigIntAlg>(a: &Mat4x4<T>, b: &Mat4x4<T>) -> bool {
        // TODO: make that constant time
        for i in 0..4 {
            for j in 0..4 {
                if a[i][j] != b[i][j] {
                    return false;
                }
            }
        }
        true
    }
}
