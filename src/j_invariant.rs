use crate::bigint::Power;

use super::polynomial;
use super::eisenstein;
use super::delta;
use super::term_builder;

/// j_invariant * q
pub fn j_invariant1(order: i64) -> polynomial::Polynomial {
    let di = delta::delta1_inverse(order);
    let e4 = eisenstein::eisenstein4(order).power(3);
    let e4 = e4.omit_high_order_q(order as i64);
    let j = di * e4;
    j.omit_high_order_q(order as i64)
}

// j_invariant
pub fn j_invariant(order: i64) -> polynomial::Polynomial {
    let j1 = j_invariant1(order + 1);
    let t = term_builder::TermBuilder::new()
            .qpow(-1)
            .build()
            .to_pol();
    j1 * t
}

#[test]
fn j_invariant_test1() {
    assert_eq_str!(j_invariant1(2), "196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(3), "21493760 q^3 + 196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(4), "864299970 q^4 + 21493760 q^3 + 196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(5), "20245856256 q^5 + 864299970 q^4 + 21493760 q^3 + 196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(6), "333202640600 q^6 + 20245856256 q^5 + 864299970 q^4 + 21493760 q^3 + 196884 q^2 + 744 q + 1");
}

#[test]
fn j_invariant_test() {
    assert_eq_str!(j_invariant(2), "21493760 q^2 + 196884 q + 744 + q^-1");
    assert_eq_str!(j_invariant(3), "864299970 q^3 + 21493760 q^2 + 196884 q + 744 + q^-1");
    assert_eq_str!(j_invariant(4), "20245856256 q^4 + 864299970 q^3 + 21493760 q^2 + 196884 q + 744 + q^-1");
    assert_eq_str!(j_invariant(5), "333202640600 q^5 + 20245856256 q^4 + 864299970 q^3 + 21493760 q^2 + 196884 q + 744 + q^-1");
    assert_eq_str!(j_invariant(6), "4252023300096 q^6 + 333202640600 q^5 + 20245856256 q^4 + 864299970 q^3 + 21493760 q^2 + 196884 q + 744 + q^-1");
}

