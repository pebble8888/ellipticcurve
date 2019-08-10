extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_integer::Integer;
use num_bigint::BigInt;
use super::polynomial;
use super::unitbuilder;
use super::division_polynomial;

use crate::bigint::{Power};

type UnitBuilder = unitbuilder::UnitBuilder;
type Polynomial = polynomial::Polynomial;

pub struct SchoofResult {
    pub l: BigInt,
    pub a: BigInt,
}

pub fn schoof(a: &BigInt, b: &BigInt, q: &BigInt) -> Vec<SchoofResult> {
    let mut schoof_result: Vec<SchoofResult> = Vec::new();
    for l in &vec![BigInt::from(3), BigInt::from(5)] {
        println!("{} l:{}", line!(), l);
        let ql = q.mod_floor(l);
        if l >= q {
            break;
        }
        println!("q:{} l:{} ql:{}", q, l, ql);
        let jmax: BigInt = (l-1) / 2;
        println!("j:[1, {}]", jmax);

        // (b) x'
        let nn1 = division_polynomial::omega(a, b, &ql);
        //println!("nn1:{} {}", line!(), nn1);

        let nn1 = nn1.reduction_modular(a, b, q);
        println!("nn1:{} {}", line!(), nn1);

        let nn2 = UnitBuilder::new().coef_i(1).ypow(&q.power_i(2)).finalize()
                * division_polynomial::psi(a, b, &ql).power_i(3);
        //println!("nn2:{} {}", line!(), nn2);

        let nn2 = nn2.reduction_modular(a, b, q);
        //println!("nn2:{} {}", line!(), nn2);

        let n1 = (nn1 - nn2).power_i(2);
        //println!("n1:{} {}", line!(), n1);

        let n1 = n1.modular(q);
        //println!("n1:{} {}", line!(), n1);

        let n2 = - (division_polynomial::phi(a, b, &ql)
                     + UnitBuilder::new().coef_i(1).xpow(&q.power_i(2)).finalize()
                       * division_polynomial::psi(a, b, &ql).power_i(2));
        //println!("n2:{} {}", line!(), n2);

        let n3 = (division_polynomial::phi(a, b, &ql)
                     - UnitBuilder::new().coef_i(1).xpow(&q.power_i(2)).finalize()
                       * division_polynomial::psi(a, b, &ql).power_i(2)).power_i(2);
        println!("n3:{} {}", line!(), n3);

        let num1 = (n1 + n2 * n3).modular(q); 
        println!("num1:{} {}", line!(), num1);

        let num1 = num1.reduction_modular(a, b, q);
        
        println!("num1:{} {}", line!(), num1);

        let den1 = (division_polynomial::psi(a, b, &ql).power_i(2))
                 * ((division_polynomial::phi(a, b, &ql) - 
                        UnitBuilder::new().coef_i(1).xpow(&q.power_i(2)).finalize()
                   * division_polynomial::psi(a, b, &ql).power_i(2)).power_i(2));
        //println!("den1:{} {}", line!(), den1);
        
        let den1 = den1.reduction_modular(a, b, q);
        println!("den1:{} {}", line!(), den1);

        let mut found = false;
        let mut jj: BigInt = BigInt::from(0);
        for j in num_iter::range(BigInt::from(1), jmax + BigInt::from(1)) {
            println!("{} j:{}", line!(), j);
            // x
            let num2 = division_polynomial::phi(a, b, &j).to_frob(q);
            //println!("num2:{} {}", line!(), num2);

            let num2 = num2.reduction_modular(a, b, q);
            println!("num2:{} {}", line!(), num2);

            let den2 = division_polynomial::psi(a, b, &j).to_frob(q).power_i(2);
            println!("den2:{} {}", line!(), den2);

            let p1 = &num1 * &den2 - &num2 * &den1;
            let p1 = p1.modular(q);
            //println!("p1:{} {}", line!(), p1);

            let p1 = p1.reduction_modular(a, b, q);
            //println!("p1:{} {}", line!(), p1);

            let p1 = p1.modular(q);
            //println!("p1:{} {}", line!(), p1);

            let psil = division_polynomial::psi(a, b, &l).modular(q);
            println!("{} psi({}):{}", line!(), l, psil);

            let p1 = p1.polynomial_modular(&psil, q);

            let p1 = &p1.modular(q);
            println!("{} pol % psi({}):{}", line!(), l, p1);

            if p1.is_zero() {
                // to iii.
                found = true;
                jj = j.clone();
                break;
            }
        }
        if found {
            // iii
            println!("x found");
            // y
            let g = UnitBuilder::new().coef_i(1).ypow(&q.power_i(2)).finalize().to_pol();
            println!("g:{} {}", line!(), g);
            let g = g.reduction_modular(a, b, q);
            println!("g:{} {}", line!(), g);

            let omg = division_polynomial::omega(a, b, &ql);
            let omg = omg.modular(q);
            //println!("omg:{} {}", line!(), omg);

            let d = omg - &g * division_polynomial::psi(a, b, &ql).power_i(3);
            //println!("d:{} {}", line!(), d);
            let d = d.reduction_modular(a, b, q);
            println!("d:{} {}", line!(), d);

            let e = - division_polynomial::phi(a, b, &ql) +
                    UnitBuilder::new().coef_i(2).xpow(&q.power_i(2)).finalize()
                    * division_polynomial::psi(a, b, &ql).power_i(2);
            let e = e.modular(q);
            //println!("e:{} {}", line!(), e);

            let e = e.reduction_modular(a, b, q);
            println!("e:{} {}", line!(), e);

            let f = division_polynomial::phi(a, b, &ql) - 
                UnitBuilder::new().coef_i(1).xpow(&q.power_i(2)).finalize() * division_polynomial::psi(a, b, &ql).power_i(2);
            let f = f.modular(q);
            //println!("f:{} {}", line!(), f);

            let f = f.reduction_modular(a, b, q);
            println!("f:{} {}", line!(), f);

            let num3 = &d * (e * f.power_i(2) - d.clone().power_i(2)) - &g * division_polynomial::psi(a, b, &ql).power_i(3) * f.power_i(3);
            //println!("num3:{} {}", line!(), num3);
            let num3 = num3.reduction_modular(a, b, q);
            println!("num3:{} {}", line!(), num3);

            let num3 = num3 * division_polynomial::psi(a, b, &ql);
            let num3 = num3.reduction_modular(a, b, q);
            println!("num3:{} {}", line!(), num3);

            let den3 = division_polynomial::psi(a, b, &ql).power_i(3) * f.power_i(3) * division_polynomial::psi(a, b, &ql);
            let den3 = den3.reduction_modular(a, b, q);
            println!("den3:{} {}", line!(), den3);

            let num4 = division_polynomial::omega(a, b, &jj).to_frob(q);
            let num4 = num4.reduction_modular(a, b, q);
            println!("num4:{} {}", line!(), num4);

            let den4 = division_polynomial::psi(a, b, &jj).to_frob(q).power_i(3);
            println!("den4:{} {}", line!(), den4);

            let den4 = den4.reduction_modular(a, b, q);
            println!("den4:{} {}", line!(), den4);

            let p7 = num3 * den4 - num4 * den3;
            //println!("p7:{} {}", line!(), p7);

            let p8 = p7.reduction_modular(a, b, q);
            println!("p8:{} {}", line!(), p8);

            let p9: Polynomial;
            if p8.has_y() {
                p9 = p8 / UnitBuilder::new().coef_i(1).ypow_i(1).finalize();
            } else {
                p9 = p8.clone();
            }
            println!("{} p9:{}", line!(), p9);

            let p9 = p9.modular(q);
            let psil = division_polynomial::psi(a, b, &l).modular(q);
            println!("{} psi({}):{}", line!(), l, psil);

            let p11 = p9.polynomial_modular(&psil, q);
            println!("{} pol % psi({}):{}", line!(), l, p11);

            if p11.is_zero() {
                schoof_result.push(SchoofResult { l: l.clone(), a: jj.clone() });
                println!("(iii) y  a = {} mod {}", jj.clone(), l);
            } else {
                schoof_result.push(SchoofResult { l: l.clone(), a: - jj.clone() });
                println!("(iii) y  a = {} mod {}", (-jj.clone()), l);
            }
        } else {
            // (d)
            println!("x not found");
            let mut found = false;
            let mut w = BigInt::from(1);
            for i in num_iter::range(BigInt::from(1), l.clone()) {
                if i.power_i(2).mod_floor(l) == ql {
                    found = true;
                    w = i.clone();
                    break;
                }
            }
            if !found {
                schoof_result.push(SchoofResult { l: l.clone(), a: BigInt::from(0) });
                println!("(d)  a = 0 mod l because of w not found (d)");
            } else {
                // (e) x
                let p12 = UnitBuilder::new().coef_i(1).xpow(&q.clone()).finalize()
                       * division_polynomial::psi(a, b, &w).power_i(2) - division_polynomial::phi(a, b, &w);
                let p12 = p12.modular(q);
                let p13 = p12.reduction_modular(a, b, q);
                println!("{} p13:{}", line!(), p13);
                let p14 = division_polynomial::psi(q, b, &l).modular(q);
                println!("{} p14:{}", line!(), p14);

                let p15 = p13.polynomial_modular(&p14, q);
                println!("{} p15 {}", line!(), p15);
                if !p15.is_zero() {
                    schoof_result.push(SchoofResult { l: l.clone(), a: BigInt::from(0) });
                    println!("(e) y  a = 0 mod {} because of gcd = 1 (e)", l);
                } else {
                    // (e) y
                    let p16 = (UnitBuilder::new().coef_i(1).ypow(&a.clone()).finalize()
                              * division_polynomial::psi(a, b, &w).power_i(3) - division_polynomial::omega(a, b, &w)) /
                            UnitBuilder::new().coef_i(1).ypow_i(1).finalize();
                    let p16 = p16.modular(q);
                    let p17 = p16.reduction_modular(a, b, q);
                    println!("{}", p17);
                    if p17.is_gcd_one(&division_polynomial::psi(a, b, &l), q) {
                        let jj = - BigInt::from(2) * w;
                        schoof_result.push(SchoofResult { l: l.clone(), a: jj.clone() });
                        println!("(e) y  a = {} mod {}", &jj, l);
                    } else {
                        let jj = BigInt::from(2) * w;
                        schoof_result.push(SchoofResult { l: l.clone(), a: jj.clone() });
                        println!("(e) y  a = {} mod {}", &jj, l);
                    }
                }
            }
        }
    }
    schoof_result
}

#[test]
fn schoof_test7() {
    let q = 7;
    let schoof_result = schoof(&BigInt::from(2), &BigInt::from(1), &BigInt::from(q));
    assert_eq!(schoof_result.len(), 2);
    assert_eq!(schoof_result[0].l, BigInt::from(3));
    assert_eq!(schoof_result[0].a, BigInt::from(0));
    assert_eq!(schoof_result[1].l, BigInt::from(5));
    assert_eq!(schoof_result[1].a, BigInt::from(-2));
}

#[test]
fn schoof_test19() {
    let q = 19;
    let schoof_result = schoof(&BigInt::from(2), &BigInt::from(1), &BigInt::from(q));
    assert_eq!(schoof_result.len(), 2);
    assert_eq!(schoof_result[0].l, BigInt::from(3));
    assert_eq!(schoof_result[0].a, BigInt::from(-1));
    assert_eq!(schoof_result[1].l, BigInt::from(5));
    assert_eq!(schoof_result[1].a, BigInt::from(-2));
}

