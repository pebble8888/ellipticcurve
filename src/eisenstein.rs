extern crate num_iter;

use num_bigint::BigInt;
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;
use super::divisor;

/// E_4(q)
pub fn eisenstein4(max_q_order: &BigInt) -> polynomial::Polynomial {
    let mut pol = term_builder::TermBuilder::new()
        .coef(&BigInt::from(1))
        .build()
        .to_pol();
    let order_plus_1 = max_q_order.clone() + BigInt::from(1);
    for n in num_iter::range(BigInt::from(1), order_plus_1) {
        let sigma = divisor::sigma_divisor(&n, &(BigInt::from(3)));
        let t = term_builder::TermBuilder::new()
            .coef(&(BigInt::from(240) * sigma))
            .xpow(&n)
            .build()
            .to_pol();
        pol += t;
    }
    pol
}

#[test]
fn eisenstein4_test1() {
    use crate::bigint::{Power};

    assert_eq_str!(eisenstein4(&BigInt::from(3)), "6720 x^3 + 2160 x^2 + 240 x + 1");
    assert_eq_str!(eisenstein4(&BigInt::from(4)), "17520 x^4 + 6720 x^3 + 2160 x^2 + 240 x + 1");
    assert_eq_str!(eisenstein4(&BigInt::from(5)), "30240 x^5 + 17520 x^4 + 6720 x^3 + 2160 x^2 + 240 x + 1");
    assert_eq_str!(eisenstein4(&BigInt::from(6)), "60480 x^6 + 30240 x^5 + 17520 x^4 + 6720 x^3 + 2160 x^2 + 240 x + 1");

    let t3 = eisenstein4(&BigInt::from(2));
    assert_eq_str!(t3.power(3), "10077696000 x^6 + 3359232000 x^5 + 387244800 x^4 + 16934400 x^3 + 179280 x^2 + 720 x + 1");

    let t4 = eisenstein4(&BigInt::from(3));
    assert_eq_str!(t4.power(3), "303464448000 x^9 + 292626432000 x^8 + 126572544000 x^7 + 31115059200 x^6 + 4607539200 x^5 + 396921600 x^4 + 16954560 x^3 + 179280 x^2 + 720 x + 1");

    let t5 = eisenstein4(&BigInt::from(4));
    assert_eq_str!(t5.power(3), "5377771008000 x^12 + 6188120064000 x^11 + 4362564096000 x^10 + 2050306560000 x^9 + 708308755200 x^8 + 181773158400 x^7 + 34369574400 x^6 + 4632768000 x^5 + 396974160 x^4 + 16954560 x^3 + 179280 x^2 + 720 x + 1");
}

