//use num_integer::Integer;
//use std::fmt;
//use super::polynomial;
//use super::term_builder::TermBuildable;
//use super::term_builder;
//use crate::bigint::Inverse;
//use num_traits::Zero;
//use num_traits::ToPrimitive;

use crate::bigint::Power;
use num_bigint::BigInt;
use num_traits::One;
use super::elliptic_curve;

#[derive(Debug, Clone)]
pub struct Secp256k1 {
    pub ec: elliptic_curve::EllipticCurve,
    pub g: elliptic_curve::ECPoint,
}

impl Secp256k1 {
    pub fn new() -> Self {
        let p = BigInt::from(2).power(256) - BigInt::from(2).power(32) - BigInt::from(977);
        let a = BigInt::from(0);
        let b = BigInt::from(7);
        let ec = elliptic_curve::EllipticCurve::new_raw(&a, &b, &p);
        let gx = BigInt::from(5) * BigInt::from(10).power((12 * 3 + 2) * 2)
                   + BigInt::from(50_662_630_222_773_436_695_787_188_951_685_343_262u128) * BigInt::from(10).power(12 * 3 + 2)
                   + BigInt::from(50_603_453_777_594_175_500_187_360_389_116_729_240u128);
        let gy = BigInt::from(3) * BigInt::from(10).power((12 * 3 + 2) * 2)
                   + BigInt::from(26_705_100_207_588_169_780_830_851_305_070_431_844u128) * BigInt::from(10).power(12 * 3 + 2)
                   + BigInt::from(71_273_380_659_243_275_938_904_335_757_337_482_424u128);
        let g = elliptic_curve::ECPoint::new(
            &gx,
            &gy,
            &One::one());

        Secp256k1 {
            ec: ec,
            g: g, 
        }
    }
}

#[test]
fn secp256k1_test1() {
    let curve = Secp256k1::new();
    let ec = curve.ec;
    let g = curve.g;

    assert!(ec.is_on_curve(&g));

    assert_eq_str!(ec.multiply_scalar(&g, &1.into()), 
        "(55066263022277343669578718895168534326250603453777594175500187360389116729240, 32670510020758816978083085130507043184471273380659243275938904335757337482424)");
    assert_eq_str!(ec.multiply_scalar(&g, &2.into()),
        "(89565891926547004231252920425935692360644145829622209833684329913297188986597, 12158399299693830322967808612713398636155367887041628176798871954788371653930)");
    assert_eq_str!(ec.multiply_scalar(&g, &3.into()),
        "(112711660439710606056748659173929673102114977341539408544630613555209775888121, 25583027980570883691656905877401976406448868254816295069919888960541586679410)");
}

#[test]
#[ignore]
fn secp256k1_test2() {
    use rand::Rng;
    use std::io::{self, Write};

    let mut rng = rand::thread_rng();

    let curve = Secp256k1::new();
    let ec = curve.ec;
    let g = curve.g;

    let sx = BigInt::from(6) * BigInt::from(10).power(38 * 2)
                        + BigInt::from(81671193674923859464776141273914116991u128) * BigInt::from(10).power(38)
                        + BigInt::from(87587739043560034303210379772692937810u128);
    let sy = BigInt::from(9) * BigInt::from(10).power(38 * 2)
                        + BigInt::from(88118691354782330710178928771578899669u128) * BigInt::from(10).power(38)
                        + BigInt::from(23567543651789733849964475999540762862u128);
    let s = elliptic_curve::ECPoint::new(
        &sx,
        &sy,
        &One::one());
    
    loop {
        let i0: u128 = rng.gen();
        let i1: u128 = rng.gen();
        let sk: BigInt = BigInt::from(i0)
            + BigInt::from(i1) * BigInt::from(2).power(128);
        //print!("sk:{}\n", sk);
        //io::stdout().flush().unwrap();
        let pk = ec.multiply_scalar(&g, &sk);    
        //print!("pk:{}\n", pk);
        print!(".");
        io::stdout().flush().unwrap();
        if pk == s {
            print!("found sk:{}", sk);
            break;
        }
    }
}

