extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;
use num_traits::{Zero};
use super::term_builder;
use super::term_builder::TermBuildable;
use super::modular_polynomial;
use super::elliptic_curve;
use super::polynomial;

type TermBuilder = term_builder::TermBuilder;

pub struct SEAResult {
    pub gcd: polynomial::Polynomial,
    pub degree_of_gcd: BigInt, 
    pub is_elkies_prime: bool,
    pub isogeny_j_invariants: Vec<BigInt>
}

/// SEA algorithm
/// TODO: in implementation
pub fn sea(ec: &elliptic_curve::EllipticCurve, l: u64) -> SEAResult {
    let mut mpol = modular_polynomial::modular_polynomial(l);
    mpol.modular_assign(&ec.p);
    let mut mpol = mpol.eval_y(&ec.j_invariant());
    mpol.modular_assign(&ec.p);
    let pol = TermBuilder::new().xpow(&ec.p.clone()).build()
            - TermBuilder::new().xpow(1).build();
    let gcd = pol.gcd(&mpol, &ec.p);
    // elkies prime for degree 1, 2, l+1
    // atkins prime for degree 0
    let degree = gcd.clone().degree_x();
    let is_elkies_prime = gcd.degree_x() > Zero::zero();
    let mut isogeny_j_invariants: Vec<BigInt> = Vec::new();
    if is_elkies_prime {
        for j in num_iter::range(BigInt::from(0), ec.p.clone()) {
            let mut tmp = gcd.clone().eval_x(&j.clone());
            tmp.modular_assign(&ec.p);
            if tmp.is_zero() {
                isogeny_j_invariants.push(j);
                if degree.clone() == BigInt::from(isogeny_j_invariants.len()) ||
                   isogeny_j_invariants.len() == 2
                {
                    break;
                }
            }
        }
    }
    SEAResult {
        gcd: gcd.clone(),
        degree_of_gcd: degree.clone(),
        is_elkies_prime: is_elkies_prime,
        isogeny_j_invariants: isogeny_j_invariants,
    }
}

#[test]
fn sea_test1() {
    let p = BigInt::from(23);
    let ec = elliptic_curve::EllipticCurve::new(&BigInt::from(1), &BigInt::from(7), &p);
    let l = 3;
    let result = sea(&ec, l);
    assert_eq_str!(result.gcd, "x^2 + 4 x + 3");

    assert_eq!(result.gcd.degree_x(), BigInt::from(2));
    assert!(result.is_elkies_prime);
    assert_eq!(result.isogeny_j_invariants.len(), 2);
    assert_eq_str!(result.isogeny_j_invariants[0], "20");
    assert_eq_str!(result.isogeny_j_invariants[1], "22");
}

#[test]
fn sea_test2() {
    let p = BigInt::from(131);
    let ec = elliptic_curve::EllipticCurve::new(&BigInt::from(1), &BigInt::from(23), &p);
    let l = 5;
    let result = sea(&ec, l);
    assert_eq_str!(result.gcd, "x^2 + 88 x + 49");
    assert_eq!(result.gcd.degree_x(), BigInt::from(2));
    assert!(result.is_elkies_prime);
    assert_eq!(result.isogeny_j_invariants.len(), 2);
    assert_eq_str!(result.isogeny_j_invariants[0], "17");
    assert_eq_str!(result.isogeny_j_invariants[1], "26");
}

