extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_integer::Integer;
use num_bigint::BigInt;
use super::polynomial;
use super::unitbuilder;
use super::unitbuilder::UnitBuildable;
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

        let nn1 = nn1.reduction_modular(a, b, q);

        let nn2 = UnitBuilder::new().coef(1).ypow(&q.power(2)).build()
                * division_polynomial::psi(a, b, &ql).power(3);

        let nn2 = nn2.reduction_modular(a, b, q);

        let mut n1 = (nn1 - nn2).power(2);

        n1.modular_assign(q);

        let n2 = - (division_polynomial::phi(a, b, &ql)
                     + UnitBuilder::new().coef(1).xpow(&q.power(2)).build()
                       * division_polynomial::psi(a, b, &ql).power(2));

        let n3 = (division_polynomial::phi(a, b, &ql)
                     - UnitBuilder::new().coef(1).xpow(&q.power(2)).build()
                       * division_polynomial::psi(a, b, &ql).power(2)).power(2);

        let mut num1 = n1 + n2 * n3;
        num1.modular_assign(q);
        println!("num1:{} {}", line!(), num1);

        let num1 = num1.reduction_modular(a, b, q);
        
        println!("num1:{} {}", line!(), num1);

        let den1 = (division_polynomial::psi(a, b, &ql).power(2))
                 * ((division_polynomial::phi(a, b, &ql) - 
                        UnitBuilder::new().coef(1).xpow(&q.power(2)).build()
                   * division_polynomial::psi(a, b, &ql).power(2)).power(2));
        
        let den1 = den1.reduction_modular(a, b, q);
        println!("den1:{} {}", line!(), den1);

        let psi_l = division_polynomial::psi(a, b, &l).modular(q);
        println!("{} psi({}):{}", line!(), l, psi_l);

        let mut found = false;
        let mut jj: BigInt = BigInt::from(0);
        for j in num_iter::range(BigInt::from(1), jmax + BigInt::from(1)) {
            println!("{} j:{}", line!(), j);
            // x
            let num2 = division_polynomial::phi(a, b, &j).to_frob(q);

            let num2 = num2.reduction_modular(a, b, q);
            println!("num2:{} {}", line!(), num2);

            let den2 = division_polynomial::psi(a, b, &j).to_frob(q).power(2);
            println!("den2:{} {}", line!(), den2);

            let mut p1 = &num1 * &den2 - &num2 * &den1;
            p1.modular_assign(q);

            let p1 = p1.reduction_modular(a, b, q);

            let mut p1 = p1.polynomial_modular(&psi_l, q);
            p1.modular_assign(q);
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
            let g = UnitBuilder::new().coef(1).ypow(&q.power(2)).build().to_pol();
            println!("g:{} {}", line!(), g);
            let g = g.reduction_modular(a, b, q);
            println!("g:{} {}", line!(), g);

            let mut omg = division_polynomial::omega(a, b, &ql);
            omg.modular_assign(q);
            //println!("omg:{} {}", line!(), omg);

            let d = omg - &g * division_polynomial::psi(a, b, &ql).power(3);
            //println!("d:{} {}", line!(), d);
            let d = d.reduction_modular(a, b, q);
            println!("d:{} {}", line!(), d);

            let mut e = - division_polynomial::phi(a, b, &ql) +
                    UnitBuilder::new().coef(2).xpow(&q.power(2)).build()
                    * division_polynomial::psi(a, b, &ql).power(2);
            e.modular_assign(q);

            let e = e.reduction_modular(a, b, q);
            println!("e:{} {}", line!(), e);

            let mut f = division_polynomial::phi(a, b, &ql) - 
                UnitBuilder::new().coef(1).xpow(&q.power(2)).build() * division_polynomial::psi(a, b, &ql).power(2);
            f.modular_assign(q);

            let f = f.reduction_modular(a, b, q);
            println!("f:{} {}", line!(), f);

            let num3 = &d * (e * f.power(2) - d.clone().power(2)) - &g * division_polynomial::psi(a, b, &ql).power(3) * f.power(3);
            let mut num3 = num3.reduction_modular(a, b, q);
            println!("num3:{} {}", line!(), num3);

            num3 *= division_polynomial::psi(a, b, &ql);
            let num3 = num3.reduction_modular(a, b, q);
            println!("num3:{} {}", line!(), num3);

            let den3 = division_polynomial::psi(a, b, &ql).power(3) * f.power(3) * division_polynomial::psi(a, b, &ql);
            let den3 = den3.reduction_modular(a, b, q);
            println!("den3:{} {}", line!(), den3);

            let num4 = division_polynomial::omega(a, b, &jj).to_frob(q);
            let num4 = num4.reduction_modular(a, b, q);
            println!("num4:{} {}", line!(), num4);

            let den4 = division_polynomial::psi(a, b, &jj).to_frob(q).power(3);

            let den4 = den4.reduction_modular(a, b, q);
            println!("den4:{} {}", line!(), den4);

            let mut p7 = num3 * den4 - num4 * den3;
            p7.modular_assign(q);

            let p8 = p7.reduction_modular(a, b, q);
            println!("p8:{} {}", line!(), p8);

            let mut p9: Polynomial;
            if p8.has_y() {
                p9 = p8 / UnitBuilder::new().coef(1).ypow(1).build();
            } else {
                p9 = p8.clone();
            }
            println!("{} p9:{}", line!(), p9);

            p9.modular_assign(q);

            let p9 = p9.polynomial_modular(&psi_l, q);
            println!("{} pol % psi({}):{}", line!(), l, p9);

            if p9.is_zero() {
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
                if i.power(2).mod_floor(l) == ql {
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
                let mut p12 = UnitBuilder::new().coef(1).xpow(&q.clone()).build()
                       * division_polynomial::psi(a, b, &w).power(2) - division_polynomial::phi(a, b, &w);
                p12.modular_assign(q);
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
                    let mut p16 = (UnitBuilder::new().coef(1).ypow(&a.clone()).build()
                              * division_polynomial::psi(a, b, &w).power(3) - division_polynomial::omega(a, b, &w)) /
                            UnitBuilder::new().coef(1).ypow(1).build();
                    p16.modular_assign(q);
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

