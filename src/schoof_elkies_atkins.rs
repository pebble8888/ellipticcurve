extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;
use super::term_builder;
use super::term_builder::TermBuildable;
use super::modular_polynomial;
use super::elliptic_curve;

type TermBuilder = term_builder::TermBuilder;

#[test]
fn sea_test1() {
    let p = BigInt::from(23);
    let ec = elliptic_curve::EllipticCurve::new(&BigInt::from(1), &BigInt::from(7), &p);
    let mut mpol = modular_polynomial::modular_polynomial(&BigInt::from(3));
    mpol.modular_assign(&p);
    let mut mpol = mpol.eval_y(&ec.j_invariant());
    mpol.modular_assign(&p);

    let pol = TermBuilder::new().coef(1).xpow(&p.clone()).build()
            - TermBuilder::new().coef(1).xpow(1).build();
    let gcd = pol.gcd(&mpol, &p);
    // elkies prime for degree 1, 2, l+1
    // atkins prime for degree 0
    assert_eq_str!(gcd, "x^2 + 4 x + 3");

    assert_eq!(gcd.degree(), BigInt::from(2));
}

#[test]
fn sea_test2() {
    let p = BigInt::from(131);
    let ec = elliptic_curve::EllipticCurve::new(&BigInt::from(1), &BigInt::from(23), &p);
    let mut mpol = modular_polynomial::modular_polynomial(&BigInt::from(5));
    mpol.modular_assign(&p);
    let mut mpol = mpol.eval_y(&ec.j_invariant());
    mpol.modular_assign(&p);
    assert_eq_str!(mpol, "x^6 + x^5 + 67 x^4 + 106 x^3 + 16 x^2 + 33 x + 41");

    let pol = TermBuilder::new().coef(1).xpow(&p.clone()).build()
            - TermBuilder::new().coef(1).xpow(1).build();
    let gcd = pol.gcd(&mpol, &p);
    assert_eq_str!(gcd, "x^2 + 88 x + 49");

    assert_eq!(gcd.degree(), BigInt::from(2));
}

