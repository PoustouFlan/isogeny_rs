// ===================================================================
// Compting general (2,2)-isogenies between theta structures
//
// NOTE: For the two steps before a product structure is reached, we
// need additional symplectic transforms which is controlled by the
// `hadamard` array of `bool`s. The purpose of these is to avoid null
// points (or dual null points) which have zero elements, which are
// incompatible with the doubling formula.
// ===================================================================

use fp2::traits::Fp2 as FqTrait;

use super::theta_point::ThetaPoint;
use super::theta_structure::ThetaStructure;
use crate::theta::theta_util::to_hadamard;

impl<Fq: FqTrait> ThetaStructure<Fq> {
    /// Returns `0xFF..FF` when <T1, T2> are isotropic and `0x00..00` otherwise.
    fn check_isotropic(T1: &ThetaPoint<Fq>, T2: &ThetaPoint<Fq>, I: &ThetaPoint<Fq>) -> u32 {
        let mut ok = u32::MAX;
        let (x1, y1, z1, t1) = T1.coords();
        let (x2, y2, z2, t2) = T2.coords();
        let (A_inv, B_inv, C_inv, D_inv) = I.coords();
        ok &= (x1 * A_inv).equals(&(y1 * B_inv));
        ok &= (z1 * C_inv).equals(&(t1 * D_inv));
        ok &= (x2 * A_inv).equals(&(z2 * C_inv));
        ok &= (y2 * B_inv).equals(&(t2 * D_inv));
        ok
    }

    /// Given the 8-torsion above the kernel, compute the codomain of the
    /// (2,2)-isogeny and the image of all points in `image_points`.
    /// Cost:
    ///   Codomain: 8S + 9M
    ///   Image: 4S + 4M per point
    pub fn two_isogeny(
        T1: &ThetaPoint<Fq>,
        T2: &ThetaPoint<Fq>,
        image_points: &mut [ThetaPoint<Fq>],
        hadamard: [bool; 2],
        verify_codomain: bool,
    ) -> (Self, u32) {
        let mut ok = u32::MAX;

        let (xA, xB, yC, yD) = T1.cond_hadamard_square_hadamard(hadamard[0]);
        let (zA, tB, zC, tD) = T2.cond_hadamard_square_hadamard(hadamard[0]);

        let xAtB = xA * tB;
        let zAxB = zA * xB;
        let zCtD = zC * tD;

        let mut A = zA * xAtB;
        let mut B = tB * zAxB;
        let mut C = zC * xAtB;
        let mut D = tD * zAxB;

        let A_inv = xB * zCtD;
        let B_inv = xA * zCtD;
        let C_inv = D;
        let D_inv = C;

        if verify_codomain {
            ok &= !T1.has_zero_coordinate();
            ok &= !T2.has_zero_coordinate();
            ok &= Self::check_isotropic(
                &ThetaPoint::from((xA, xB, yC, yD)),
                &ThetaPoint::from((zA, tB, zC, tD)),
                &ThetaPoint::from((A_inv, B_inv, C_inv, D_inv)),
            );
        }

        if hadamard[1] {
            (A, B, C, D) = to_hadamard(&A, &B, &C, &D);
        }

        let codomain = Self::new_from_coords(&A, &B, &C, &D);

        for P in image_points.iter_mut() {
            let (mut XX, mut YY, mut ZZ, mut TT) = P.cond_hadamard_square_hadamard(hadamard[0]);
            XX *= A_inv;
            YY *= B_inv;
            ZZ *= C_inv;
            TT *= D_inv;
            if hadamard[1] {
                (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
            }
            *P = ThetaPoint::from((XX, YY, ZZ, TT));
        }

        (codomain, ok)
    }

    /// Given the 4-torsion above the kernel, compute the codomain of the
    /// (2,2)-isogeny and the image of all points in `image_points`.
    /// Cost:
    ///   Codomain: 8S + 17M + 2Sqrt
    pub fn two_isogeny_4_torsion(
        &self,
        T1_prime: &ThetaPoint<Fq>,
        image_points: &mut [ThetaPoint<Fq>],
        hadamard: [bool; 2],
        verify_codomain: bool,
    ) -> (Self, u32) {
        let mut ok = u32::MAX;

        // Extract squared coordinates
        let (xAB, _, xCD, _) = T1_prime.cond_hadamard_square_hadamard(hadamard[0]);
        let (A2, B2, C2, D2) = self.null_point().cond_hadamard_square_hadamard(hadamard[0]);

        // Compute square roots required for the 4-torsion codomain structure
        let (AB, r1) = (A2 * B2).sqrt();
        let (AC, r2) = (A2 * C2).sqrt();

        ok &= r1 & r2;

        // Compute codomain coordinates
        let mut B = AB * AC;
        let mut D_inv = B * xCD;
        B *= xAB;

        let D = xCD * AB * A2;

        let mut A = xAB * A2;
        let C = A * C2;
        A *= AC;

        // Compute coordinate inverses for mapping image points
        let mut A_inv = xAB * D2;
        let mut C_inv = A_inv * B2;
        A_inv *= C2;
        let B_inv = A_inv * AB;
        A_inv *= B2;
        C_inv *= AC;
        D_inv *= B2;

        let mut codomain_A = A;
        let mut codomain_B = B;
        let mut codomain_C = C;
        let mut codomain_D = D;

        if hadamard[1] {
            (codomain_A, codomain_B, codomain_C, codomain_D) = to_hadamard(&codomain_A, &codomain_B, &codomain_C, &codomain_D);
        }

        if verify_codomain {
            ok &= !codomain_A.is_zero();
        }

        let codomain = Self::new_from_coords(&codomain_A, &codomain_B, &codomain_C, &codomain_D);

        // Push image points through the isogeny
        for P in image_points.iter_mut() {
            let (mut XX, mut YY, mut ZZ, mut TT) = P.cond_hadamard_square_hadamard(hadamard[0]);
            XX *= A_inv;
            YY *= B_inv;
            ZZ *= C_inv;
            TT *= D_inv;
            if hadamard[1] {
                (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
            }
            *P = ThetaPoint::new(&XX, &YY, &ZZ, &TT);
        }

        (codomain, ok)
    }

    /// Computes a (2,2)-isogeny from a level 2 theta structure using a
    /// 2-torsion kernel point (the null point).
    /// Cost:
    ///   Codomain: 4S + 10M + 3Sqrt
    pub fn two_isogeny_2_torsion(
        &self,
        image_points: &mut [ThetaPoint<Fq>],
        hadamard: [bool; 2],
        verify_codomain: bool,
    ) -> (Self, u32) {
        let mut ok = u32::MAX;

        // Extract squared coordinates of the null point
        let (A2, B2, C2, D2) = self.null_point().cond_hadamard_square_hadamard(hadamard[0]);

        // Compute square roots required for the 2-torsion codomain structure
        let (val_B, r1) = (A2 * B2).sqrt();
        let (val_C, r2) = (A2 * C2).sqrt();
        let (val_D, r3) = (A2 * D2).sqrt();

        ok &= r1 & r2 & r3;

        // Compute coordinate inverses for mapping image points
        let mut A_inv = C2 * D2;
        let B_inv = A_inv * val_B;
        A_inv *= B2;
        let C_inv = D2 * B2 * val_C;
        let D_inv = C2 * B2 * val_D;

        let mut codomain_A = A2;
        let mut codomain_B = val_B;
        let mut codomain_C = val_C;
        let mut codomain_D = val_D;

        if hadamard[1] {
            (codomain_A, codomain_B, codomain_C, codomain_D) = to_hadamard(&codomain_A, &codomain_B, &codomain_C, &codomain_D);
        }

        if verify_codomain {
            ok &= !codomain_A.is_zero();
        }

        let codomain = Self::new_from_coords(&codomain_A, &codomain_B, &codomain_C, &codomain_D);

        // Push image points through the isogeny
        for P in image_points.iter_mut() {
            let (mut XX, mut YY, mut ZZ, mut TT) = P.cond_hadamard_square_hadamard(hadamard[0]);
            XX *= A_inv;
            YY *= B_inv;
            ZZ *= C_inv;
            TT *= D_inv;
            if hadamard[1] {
                (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
            }
            *P = ThetaPoint::new(&XX, &YY, &ZZ, &TT);
        }

        (codomain, ok)
    }
}
