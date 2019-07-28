//extern crate num;
extern crate num_bigint; extern crate num_traits;
extern crate num_iter;

//use num::BigInt;
use num::pow::pow;
//use num_bigint::pow;
use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use super::polynomial;
use super::unit;
use super::unitbuilder;

type Unit = unit::Unit;
type UnitBuilder = unitbuilder::UnitBuilder;
type Polynomial = polynomial::Polynomial;

pub fn psi(a: &BigInt, b: &BigInt, n: &BigInt) -> Polynomial {
    assert!(*n >= Zero::zero());
    if *n == Zero::zero() {
        return Polynomial::new();
    } else if *n == BigInt::from(1) {
        return UnitBuilder::new().coef(1).finalize().to_pol();
    } else if *n == BigInt::from(2) {
        return UnitBuilder::new().coef(2).ypow(1).finalize().to_pol();
    } else if *n == BigInt::from(3) {
        Polynomial { units: vec![ 
             UnitBuilder::new().coef(3).xpow(4).finalize(),
             Unit { coef: 6 * a,  xpow: BigInt::from(2),  ypow: Zero::zero() },
             Unit { coef: 12 * b, xpow: One::one(),       ypow: Zero::zero() },
             Unit { coef: - pow(a.clone(),2), xpow: Zero::zero(), ypow: Zero::zero() }, ] }
    } else if *n == BigInt::from(4) {
        UnitBuilder::new().coef(4).ypow(1).finalize().to_pol() *
        Polynomial { units: vec![
            Unit { coef: One::one(),    xpow: BigInt::from(6), ypow: Zero::zero(), },
            Unit { coef: 5 * a,         xpow: BigInt::from(4), ypow: Zero::zero(), },
            Unit { coef: 20 * b,        xpow: BigInt::from(3), ypow: Zero::zero(), },
            Unit { coef: -5 * pow(a.clone(),2), xpow: BigInt::from(2), ypow: Zero::zero(), },
            Unit { coef: -4 * a * b,    xpow: One::one(),      ypow: Zero::zero(), },
            Unit { coef: BigInt::from(-8) * pow(b.clone(),2) - pow(a.clone(),3), xpow: Zero::zero(), ypow: Zero::zero(), },
        ] }
    } else if n % BigInt::from(2) == One::one() {
        let m: BigInt = (n-1) / 2;
        assert!(&m < n);
        let e = psi(a, b, &(m.clone()+2));
        let f = psi(a, b, &(m.clone())).power(&BigInt::from(3));
        let g = psi(a, b, &(m.clone()-1));
        let h = psi(a, b, &(m.clone()+1)).power(&BigInt::from(3));
        (e*f - g*h).ec_reduction(a.clone(), b.clone())
    } else {
        let m: BigInt = n / 2;
        assert!(&m < n);
        let e = psi(a, b, &(m.clone()+2));
        let f = psi(a, b, &(m.clone()-1)).power(&BigInt::from(2));
        let g = psi(a, b, &(m.clone()-2));
        let h = psi(a, b, &(m.clone()+1)).power(&BigInt::from(2));
        let i = psi(a, b, &m.clone()) * (e*f - g*h);
        let j = i / UnitBuilder::new().coef(2).ypow(1).finalize().to_pol();
        j.ec_reduction(a.clone(), b.clone())
    }
}

pub fn phi(a: &BigInt, b: &BigInt, n: &BigInt) -> Polynomial {
    assert!(*n >= One::one());
    let i = UnitBuilder::new().coef(1).xpow(1).finalize().to_pol() * psi(a, b, &(n.clone())).power(&BigInt::from(2))
              - psi(a, b, &(n.clone() + 1)) * psi(a, b, &(n.clone()-1));
    i.ec_reduction(a.clone(), b.clone())
}

pub fn omega(a: &BigInt, b: &BigInt, n: &BigInt) -> Polynomial {
    assert!(*n >= One::one());
    if *n == One::one() {
        UnitBuilder::new().coef(1).ypow(1).finalize().to_pol()
    } else {
        let i = (psi(a, b, &(n.clone()+2)) * psi(a, b, &(n.clone()-1)).power(&BigInt::from(2))
                     - psi(a, b, &(n.clone()-2)) * psi(a, b, &(n.clone()+1)).power(&BigInt::from(2)))
                  / UnitBuilder::new().coef(4).ypow(1).finalize().to_pol();
        i.ec_reduction(a.clone(), b.clone())
    }
}

#[test]
fn division_polynomial_test() {
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

    // phi
    let phi1 = phi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(1));
    assert_eq!(phi1.to_string(), "x");

    let phi2 = phi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(2));
    assert_eq!(phi2.to_string(), "x^4 - 2 x^2 - 8 x + 1");
    
    let phi3 = phi(&BigInt::from(1), &BigInt::from(1), &BigInt::from(3));
    assert_eq!(phi3.to_string(), "x^9 - 12 x^7 - 96 x^6 + 30 x^5 - 24 x^4 + 84 x^3 + 48 x^2 + 105 x + 72");

    // omega
    let omega1 = omega(&BigInt::from(1), &BigInt::from(1), &BigInt::from(1));
    assert_eq!(omega1.to_string(), "y");

    let omega2 = omega(&BigInt::from(1), &BigInt::from(1), &BigInt::from(2));
    assert_eq!(omega2.to_string(), "x^6 + 5 x^4 + 20 x^3 - 5 x^2 - 4 x - 9");
}

