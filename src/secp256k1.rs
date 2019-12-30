
#[test]
#[ignore]
fn secp256k1_test1() {
    use crate::bigint::Power;
    use num_bigint::BigInt;
    //use num_integer::Integer;
    //use std::fmt;
    //use std::vec;
    //use std::ops::Deref;
    //use super::polynomial;
    //use super::term_builder::TermBuildable;
    //use super::term_builder;
    //use crate::bigint::Inverse;
    //use num_traits::Zero;
    use num_traits::One;
    //use num_traits::ToPrimitive;
    use super::elliptic_curve;

    let p = BigInt::from(2).power(256) - BigInt::from(2).power(32) - BigInt::from(977);
    let a = BigInt::from(0);
    let b = BigInt::from(7);
    let ec = elliptic_curve::EllipticCurve::new_raw(&a, &b, &p);

    let base_x = BigInt::from(5) * BigInt::from(10).power((12 * 3 + 2) * 2)
               + BigInt::from(50_662_630_222_773_436_695_787_188_951_685_343_262u128) * BigInt::from(10).power(12 * 3 + 2)
               + BigInt::from(50_603_453_777_594_175_500_187_360_389_116_729_240u128);
    let base_y = BigInt::from(3) * BigInt::from(10).power((12 * 3 + 2) * 2)
               + BigInt::from(26_705_100_207_588_169_780_830_851_305_070_431_844u128) * BigInt::from(10).power(12 * 3 + 2)
               + BigInt::from(71_273_380_659_243_275_938_904_335_757_337_482_424u128);
    let base = elliptic_curve::ECPoint::new(
        &base_x,
        &base_y,
        &One::one(),
    );

    assert!(ec.is_on_curve(&base));

    assert_eq_str!(ec.multiply_scalar(&base, &1.into()), 
        "(55066263022277343669578718895168534326250603453777594175500187360389116729240, 32670510020758816978083085130507043184471273380659243275938904335757337482424)");
    assert_eq_str!(ec.multiply_scalar(&base, &2.into()),
        "");
}

