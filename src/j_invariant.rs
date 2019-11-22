use num_bigint::BigInt;
use crate::bigint::Power;
use crate::term_builder::TermBuildable;

use super::polynomial;
use super::eisenstein;
use super::delta;
use super::term_builder;

/// j_invariant * q
pub fn j_invariant1(order: &BigInt) -> polynomial::Polynomial {
    let di = delta::delta1_inverse(order);
    let e4 = eisenstein::eisenstein4(order).power(3);
    let e4 = e4.omit_high_order_q(&order);
    let j = di * e4;
    j.omit_high_order_q(order)
}

// j_invariant
pub fn j_invariant(order: &BigInt) -> polynomial::Polynomial {
    let j1 = j_invariant1(&(order + BigInt::from(1)));
    let t = term_builder::TermBuilder::new()
            .coef(&BigInt::from(1))
            .qpow(&BigInt::from(-1))
            .build()
            .to_pol();
    j1 * t
}

#[test]
fn j_invariant_test1() {
    assert_eq_str!(j_invariant1(&BigInt::from(2)), "196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(3)), "21493760 q^3 + 196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(4)), "864299970 q^4 + 21493760 q^3 + 196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(5)), "20245856256 q^5 + 864299970 q^4 + 21493760 q^3 + 196884 q^2 + 744 q + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(6)), "333202640600 q^6 + 20245856256 q^5 + 864299970 q^4 + 21493760 q^3 + 196884 q^2 + 744 q + 1");
    // only for relase build
    //assert_eq_str!(j_invariant1(&BigInt::from(30)), "13798375834642999925542288376 q^30 + 4365689224858876634610401280 q^29 + 1353563541518646878675077500 q^28 + 410789960190307909157638144 q^27 + 121883284330422510433351500 q^26 + 35307453186561427099877376 q^25 + 9971041659937182693533820 q^24 + 2740630712513624654929920 q^23 + 731811377318137519245696 q^22 + 189449976248893390028800 q^21 + 47438786801234168813250 q^20 + 11459912788444786513920 q^19 + 2662842413150775245160 q^18 + 593121772421445058560 q^17 + 126142916465781843075 q^16 + 25497827389410525184 q^15 + 4872010111798142520 q^14 + 874313719685775360 q^13 + 146211911499519294 q^12 + 22567393309593600 q^11 + 3176440229784420 q^10 + 401490886656000 q^9 + 44656994071935 q^8 + 4252023300096 q^7 + 333202640600 q^6 + 20245856256 q^5 + 864299970 q^4 + 21493760 q^3 + 196884 q^2 + 744 q + 1");
}

#[test]
fn j_invariant_test() {
    assert_eq_str!(j_invariant(&BigInt::from(2)), "21493760 q^2 + 196884 q + 744 + q^-1");
    assert_eq_str!(j_invariant(&BigInt::from(3)), "864299970 q^3 + 21493760 q^2 + 196884 q + 744 + q^-1");
}

