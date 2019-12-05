extern crate num_iter;

use num_bigint::BigInt;
use num_traits::One;
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;
use super::divisor;

/// E_4(q)
pub fn eisenstein4(max_q_order: i32) -> polynomial::Polynomial {
    let mut pol = polynomial::Polynomial::one();
    for n in num_iter::range(1, max_q_order + 1) {
        let sigma = divisor::sigma_divisor(n, 3);
        let t = term_builder::TermBuilder::new()
            .coef(BigInt::from(240) * sigma)
            .qpow(n)
            .build()
            .to_pol();
        pol += t;
    }
    pol
}

#[test]
fn eisenstein4_test1() {
    use crate::bigint::{Power};

    assert_eq_str!(eisenstein4(3), "6720 q^3 + 2160 q^2 + 240 q + 1");
    assert_eq_str!(eisenstein4(4), "17520 q^4 + 6720 q^3 + 2160 q^2 + 240 q + 1");
    assert_eq_str!(eisenstein4(5), "30240 q^5 + 17520 q^4 + 6720 q^3 + 2160 q^2 + 240 q + 1");
    assert_eq_str!(eisenstein4(6), "60480 q^6 + 30240 q^5 + 17520 q^4 + 6720 q^3 + 2160 q^2 + 240 q + 1");

    let t3 = eisenstein4(2);
    assert_eq_str!(t3.power(3), "10077696000 q^6 + 3359232000 q^5 + 387244800 q^4 + 16934400 q^3 + 179280 q^2 + 720 q + 1");

    let t4 = eisenstein4(3);
    assert_eq_str!(t4.power(3), "303464448000 q^9 + 292626432000 q^8 + 126572544000 q^7 + 31115059200 q^6 + 4607539200 q^5 + 396921600 q^4 + 16954560 q^3 + 179280 q^2 + 720 q + 1");

    let t5 = eisenstein4(4);
    assert_eq_str!(t5.power(3), "5377771008000 q^12 + 6188120064000 q^11 + 4362564096000 q^10 + 2050306560000 q^9 + 708308755200 q^8 + 181773158400 q^7 + 34369574400 q^6 + 4632768000 q^5 + 396974160 q^4 + 16954560 q^3 + 179280 q^2 + 720 q + 1");
}

