extern crate num_iter;

use num_bigint::BigInt;
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;
use crate::bigint::{Power};

/// (1-q)^24 * (1-q^2)^24 * .. * (1-q^order)^24
pub fn delta1(order: &BigInt) -> polynomial::Polynomial {
    let one = polynomial::Polynomial::one();
    let mut pol = term_builder::TermBuilder::new()
        .build()
        .to_pol();
    let order_plus_1 = order.clone() + BigInt::from(1);
    for n in num_iter::range(BigInt::from(1), order_plus_1) {
        let t = term_builder::TermBuilder::new().qpow(&n).build().to_pol();
        let u = (one.clone() - t).power(24);
        let u = u.omit_high_order_q(&order);
        pol *= u;
        pol = pol.omit_high_order_q(&order);
    }
    pol
}

// 1/d = 1 + (1-d) + (1-d)^2 + ...
pub fn delta1_inverse(order: &BigInt) -> polynomial::Polynomial {
    let a = polynomial::Polynomial::one() - delta1(&order.clone());
    let mut pol = polynomial::Polynomial::one();
    let order_plus_1 = order.clone() + BigInt::from(1);
    for n in num_iter::range(BigInt::from(1), order_plus_1) {
        let t = a.clone().power(&n);
        let t = t.omit_high_order_q(&order);
        pol += t;
    }
    pol
}

#[test]
fn delta1_test1() {
    assert_eq_str!(delta1(&BigInt::from(1)), "- 24 q + 1");
    assert_eq_str!(delta1(&BigInt::from(2)), "252 q^2 - 24 q + 1");
    assert_eq_str!(delta1(&BigInt::from(3)), "- 1472 q^3 + 252 q^2 - 24 q + 1");

    assert_eq_str!(delta1_inverse(&BigInt::from(1)), "24 q + 1");
    assert_eq_str!(delta1_inverse(&BigInt::from(2)), "324 q^2 + 24 q + 1");
    assert_eq_str!(delta1_inverse(&BigInt::from(3)), "3200 q^3 + 324 q^2 + 24 q + 1");
}

