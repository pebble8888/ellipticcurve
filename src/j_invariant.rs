use num_bigint::BigInt;
use crate::bigint::Power;

use super::polynomial;
use super::eisenstein;
use super::delta;

/// j_invariant * q
pub fn j_invariant1(order: &BigInt) -> polynomial::Polynomial {
    let di = delta::delta1_inverse(order);
    let e4 = eisenstein::eisenstein4(order).power(3);
    let j = di * e4;
    j.omit_high_order(order)
}

#[test]
fn j_invariant_test1() {
    assert_eq_str!(j_invariant1(&BigInt::from(2)), "196884 x^2 + 744 x + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(3)), "21493760 x^3 + 196884 x^2 + 744 x + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(4)), "864299970 x^4 + 21493760 x^3 + 196884 x^2 + 744 x + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(5)), "20245856256 x^5 + 864299970 x^4 + 21493760 x^3 + 196884 x^2 + 744 x + 1");
}

