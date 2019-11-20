use num_bigint::BigInt;
use crate::bigint::Power;

use super::polynomial;
use super::eisenstein;
use super::delta;

/// j_invariant * q
pub fn j_invariant1(order: &BigInt) -> polynomial::Polynomial {
    let di = delta::delta1_inverse(order);
    let e4 = eisenstein::eisenstein4(order).power(3);
    let e4 = e4.omit_high_order(&order);
    let j = di * e4;
    j.omit_high_order(order)
}

#[test]
fn j_invariant_test1() {
    assert_eq_str!(j_invariant1(&BigInt::from(2)), "196884 x^2 + 744 x + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(3)), "21493760 x^3 + 196884 x^2 + 744 x + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(4)), "864299970 x^4 + 21493760 x^3 + 196884 x^2 + 744 x + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(5)), "20245856256 x^5 + 864299970 x^4 + 21493760 x^3 + 196884 x^2 + 744 x + 1");
    assert_eq_str!(j_invariant1(&BigInt::from(6)), "333202640600 x^6 + 20245856256 x^5 + 864299970 x^4 + 21493760 x^3 + 196884 x^2 + 744 x + 1");
    // only for relase build
    //assert_eq_str!(j_invariant1(&BigInt::from(30)), "13798375834642999925542288376 x^30 + 4365689224858876634610401280 x^29 + 1353563541518646878675077500 x^28 + 410789960190307909157638144 x^27 + 121883284330422510433351500 x^26 + 35307453186561427099877376 x^25 + 9971041659937182693533820 x^24 + 2740630712513624654929920 x^23 + 731811377318137519245696 x^22 + 189449976248893390028800 x^21 + 47438786801234168813250 x^20 + 11459912788444786513920 x^19 + 2662842413150775245160 x^18 + 593121772421445058560 x^17 + 126142916465781843075 x^16 + 25497827389410525184 x^15 + 4872010111798142520 x^14 + 874313719685775360 x^13 + 146211911499519294 x^12 + 22567393309593600 x^11 + 3176440229784420 x^10 + 401490886656000 x^9 + 44656994071935 x^8 + 4252023300096 x^7 + 333202640600 x^6 + 20245856256 x^5 + 864299970 x^4 + 21493760 x^3 + 196884 x^2 + 744 x + 1");
}

