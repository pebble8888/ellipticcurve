extern crate num_iter;
use num_traits::One;
use super::polynomial;
use super::term_builder;

/// (1-q)^24 * (1-q^2)^24 * .. * (1-q^order)^24
pub fn delta1(order: i32) -> polynomial::Polynomial {
    let one = polynomial::Polynomial::one();
    let mut pol = polynomial::Polynomial::one();
    for n in num_iter::range(1, order + 1) {
        let t = term_builder::TermBuilder::new().qpow(n).build().to_pol();
        let u = (one.clone() - t).power_omit_high_order_q(24, order);
        pol *= u;
        pol.omit_high_order_q(order);
    }
    pol
}

// 1/d = 1 + (1-d) + (1-d)^2 + ...
pub fn delta1_inverse(order: i32) -> polynomial::Polynomial {
    let a = polynomial::Polynomial::one() - delta1(order);
    let mut pol = polynomial::Polynomial::one();
    for n in num_iter::range(1, order + 1) {
        let t = a.clone().power_omit_high_order_q(n, order);
        pol += t;
    }
    pol
}

#[test]
fn delta1_test1() {
    assert_eq_str!(delta1(1), "- 24 q + 1");
    assert_eq_str!(delta1(2), "252 q^2 - 24 q + 1");
    assert_eq_str!(delta1(3), "- 1472 q^3 + 252 q^2 - 24 q + 1");

    assert_eq_str!(delta1_inverse(1), "24 q + 1");
    assert_eq_str!(delta1_inverse(2), "324 q^2 + 24 q + 1");
    assert_eq_str!(delta1_inverse(3), "3200 q^3 + 324 q^2 + 24 q + 1");
}

