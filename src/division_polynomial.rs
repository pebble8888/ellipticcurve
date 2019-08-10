extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use super::polynomial;
use super::unitbuilder;
use super::bigint::{DivFloor, RemFloor, Power};

type UnitBuilder = unitbuilder::UnitBuilder;
type Polynomial = polynomial::Polynomial;

pub fn psi(a: &BigInt, b: &BigInt, n: &BigInt) -> Polynomial {
    assert!(*n >= Zero::zero());
    if *n == Zero::zero() {
        return Polynomial::new();
    } else if *n == BigInt::from(1) {
        return UnitBuilder::new().coef_i(1).finalize().to_pol();
    } else if *n == BigInt::from(2) {
        return UnitBuilder::new().coef_i(2).ypow_i(1).finalize().to_pol();
    } else if *n == BigInt::from(3) {
        UnitBuilder::new().coef_i(3).xpow_i(4).finalize()
        + UnitBuilder::new().coef(&(6 * a)).xpow_i(2).finalize()
        + UnitBuilder::new().coef(&(12 * b)).xpow_i(1).finalize()
        + UnitBuilder::new().coef(&( - a.clone().power_i(2))).finalize()
    } else if *n == BigInt::from(4) {
        UnitBuilder::new().coef_i(4).ypow_i(1).finalize() *
              ( UnitBuilder::new().coef_i(1).xpow_i(6).finalize()
              + UnitBuilder::new().coef(&(5 *a)).xpow_i(4).finalize()
              + UnitBuilder::new().coef(&(20 * b)).xpow_i(3).finalize()
              + UnitBuilder::new().coef(&(-5 * a.clone().power_i(2))).xpow_i(2).finalize()
              + UnitBuilder::new().coef(&(-4 * a * b)).xpow_i(1).finalize()
              + UnitBuilder::new().coef(&(BigInt::from(-8) * b.clone().power_i(2) - a.clone().power_i(3))).finalize()
              ) 
        
    } else if n.rem_floor(&BigInt::from(2)) == One::one() {
        let m: BigInt = (n-BigInt::from(1)).div_floor(&BigInt::from(2));
        assert!(&m < n);
        let e = psi(a, b, &(m.clone()+2));
        let f = psi(a, b, &(m.clone())).power_i(3);
        let g = psi(a, b, &(m.clone()-1));
        let h = psi(a, b, &(m.clone()+1)).power_i(3);
        (e*f - g*h).reduction(a, b)
    } else {
        let m: BigInt = n.div_floor(&BigInt::from(2));
        assert!(&m < n);
        let e = psi(a, b, &(m.clone()+2));
        let f = psi(a, b, &(m.clone()-1)).power_i(2);
        let g = psi(a, b, &(m.clone()-2));
        let h = psi(a, b, &(m.clone()+1)).power_i(2);
        let i = psi(a, b, &m.clone()) * (e*f - g*h);
        let j = i / UnitBuilder::new().coef_i(2).ypow_i(1).finalize();
        j.reduction(a, b)
    }
}

pub fn phi(a: &BigInt, b: &BigInt, n: &BigInt) -> Polynomial {
    assert!(*n >= One::one());
    let i = UnitBuilder::new().coef_i(1).xpow_i(1).finalize() * psi(a, b, &(n.clone())).power_i(2)
              - psi(a, b, &(n + 1)) * psi(a, b, &(n.clone()-1));
    i.reduction(a, b)
}

pub fn omega(a: &BigInt, b: &BigInt, n: &BigInt) -> Polynomial {
    assert!(*n >= One::one());
    if *n == One::one() {
        UnitBuilder::new().coef_i(1).ypow_i(1).finalize().to_pol()
    } else {
        let i = (psi(a, b, &(n.clone()+2)) * psi(a, b, &(n.clone()-1)).power_i(2)
                     - psi(a, b, &(n.clone()-2)) * psi(a, b, &(n.clone()+1)).power_i(2))
                  / UnitBuilder::new().coef_i(4).ypow_i(1).finalize();
        i.reduction(a, b)
    }
}

#[test]
fn division_polynomial_test_psi() {
    // psi
    let psi0 = psi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(0));
    assert_eq!(psi0.to_string(), "0");

    let psi1 = psi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(1));
    assert_eq!(psi1.to_string(), "1");

    let psi2 = psi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(2));
    assert_eq!(psi2.to_string(), "2 y");

    let psi3 = psi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(3));
    assert_eq!(psi3.to_string(), "3 x^4 + 6 x^2 + 12 x - 1");

    let psi4 = psi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(4));
    assert_eq!(psi4.to_string(), "4 x^6 y + 20 x^4 y + 80 x^3 y - 20 x^2 y - 16 x y - 36 y");

    let psi5 = psi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(5));
    assert_eq!(psi5.to_string(), "5 x^12 + 62 x^10 + 380 x^9 - 105 x^8 + 240 x^7 - 540 x^6 - 696 x^5 - 2045 x^4 - 1680 x^3 - 290 x^2 - 740 x - 287");
}

#[test]
fn division_polynomial_test_phi() {
    // phi
    let phi1 = phi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(1));
    assert_eq!(phi1.to_string(), "x");

    let phi2 = phi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(2));
    assert_eq!(phi2.to_string(), "x^4 - 2 x^2 - 8 x + 1");
    
    let phi3 = phi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(3));
    assert_eq!(phi3.to_string(), "x^9 - 12 x^7 - 96 x^6 + 30 x^5 - 24 x^4 + 84 x^3 + 48 x^2 + 105 x + 72");
}

#[test]
fn division_polynomial_test_omega() {
    // omega
    let omega1 = omega(&BigInt::from(1), &BigInt::from(1), &BigInt::from(1));
    assert_eq!(omega1.to_string(), "y");

    let omega2 = omega(&BigInt::from(1), &BigInt::from(1), &BigInt::from(2));
    assert_eq!(omega2.to_string(), "x^6 + 5 x^4 + 20 x^3 - 5 x^2 - 4 x - 9");

    let omega3 = omega(&BigInt::from(1), &BigInt::from(1), &BigInt::from(3));
    assert_eq!(omega3.to_string(), "x^12 y + 22 x^10 y + 220 x^9 y - 165 x^8 y - 528 x^7 y - 1868 x^6 y + 264 x^5 y - 1145 x^4 y - 400 x^3 y - 714 x^2 y - 1028 x y - 611 y");

    let omega4 = omega(&BigInt::from(1), &BigInt::from(1), &BigInt::from(4));
    assert_eq!(omega4.to_string(), "x^24 + 68 x^22 + 1232 x^21 - 1694 x^20 - 9856 x^19 - 61964 x^18 + 10032 x^17 - 170561 x^16 - 23040 x^15 - 341752 x^14 - 632800 x^13 - 1329636 x^12 - 2539264 x^11 - 4549752 x^10 - 4618144 x^9 - 5622449 x^8 - 6799872 x^7 - 4564812 x^6 - 4212208 x^5 - 3158750 x^4 - 892544 x^3 - 248700 x^2 - 222864 x - 40879");

    let omega5 = omega(&BigInt::from(1), &BigInt::from(1), &BigInt::from(5));
    assert_eq!(omega5.to_string(), "x^36 y + 162 x^34 y + 4692 x^33 y - 10659 x^32 y - 107712 x^31 y - 902224 x^30 y + 556512 x^29 y - 3417068 x^28 y + 2557376 x^27 y - 24924744 x^26 y - 69151824 x^25 y - 257703372 x^24 y - 686331072 x^23 y - 1968515376 x^22 y - 2825185248 x^21 y - 5467087026 x^20 y - 10374222912 x^19 y - 12843672372 x^18 y - 23464263816 x^17 y - 28086809658 x^16 y - 28443733056 x^15 y - 33582134832 x^14 y - 29309513952 x^13 y - 19226935196 x^12 y - 19770442944 x^11 y - 13785051976 x^10 y - 12217620304 x^9 y - 14642004444 x^8 y - 14782274112 x^7 y - 13037393232 x^6 y - 8387833632 x^5 y - 2784562631 x^4 y + 221827904 x^3 y + 446882082 x^2 y + 112442324 x y + 30699397 y");
}

