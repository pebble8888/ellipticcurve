extern crate num_iter;

use num_bigint::BigInt;
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;

use crate::bigint::{Power};

/// prod (1-q^n)^24
/// order: q power stable for order - 1
pub fn delta1(order: &BigInt) -> polynomial::Polynomial {
    let one = polynomial::Polynomial::one();
    let mut pol = term_builder::TermBuilder::new()
        .coef(&BigInt::from(1))
        .build()
        .to_pol();
    let order_plus_1 = order.clone() + BigInt::from(1);
    for n in num_iter::range(BigInt::from(1), order_plus_1) {
        let t = term_builder::TermBuilder::new()
            .coef(&BigInt::from(1))
            .xpow(&n)
            .build()
            .to_pol();
        let u = (one.clone() - t).power(24);
        let u = u.omit_high_order(&order);
        pol *= u;
        pol = pol.omit_high_order(&order);
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
        let t = t.omit_high_order(&order);
        pol += t;
    }
    pol
}

#[test]
fn delta1_test1() {
    assert_eq_str!(delta1(&BigInt::from(1)), "- 24 x + 1");
    assert_eq_str!(delta1(&BigInt::from(2)), "252 x^2 - 24 x + 1");
    assert_eq_str!(delta1(&BigInt::from(3)), "- 1472 x^3 + 252 x^2 - 24 x + 1");

    assert_eq_str!(delta1_inverse(&BigInt::from(1)), "24 x + 1");
    assert_eq_str!(delta1_inverse(&BigInt::from(2)), "324 x^2 + 24 x + 1");
    assert_eq_str!(delta1_inverse(&BigInt::from(3)), "3200 x^3 + 324 x^2 + 24 x + 1");
}

