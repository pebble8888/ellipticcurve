extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use super::polynomial;
use super::unit;
use super::unitbuilder;
use super::division_polynomial;

use crate::bigint::{Power, RemFloor};

type Unit = unit::Unit;
type UnitBuilder = unitbuilder::UnitBuilder;
type Polynomial = polynomial::Polynomial;

pub fn schoof(a: &BigInt, b: &BigInt, q: &BigInt) {
    for l in &vec![BigInt::from(3), BigInt::from(5)] {
        println!("{} l:{}", line!(), l.to_string());
        let ql = q.rem_floor(l);
        if l >= q {
            break;
        }
        println!("q:{} l:{} ql:{}", q.to_string(), l.to_string(), ql.to_string());
        let jmax: BigInt = (l-1) / 2;
        println!("j:[1, {}]", jmax.to_string());

        // (b) x'
        let nn1 = division_polynomial::omega(a, b, &ql);
        //println!("nn1:{} {}", line!(), nn1.to_string());

        let nn2 = Unit {
                         coef: One::one(),
                         xpow: Zero::zero(),
                         ypow: q.power(&BigInt::from(2)),
                       }.to_pol() * division_polynomial::psi(a, b, &ql).power(&BigInt::from(3));
        //println!("nn2:{} {}", line!(), nn2.to_string());

        let n1 = (nn1 - nn2).power(&BigInt::from(2));
        //println!("n1:{} {}", line!(), n1.to_string());

        let n2 = - (division_polynomial::phi(a, b, &ql)
                     + Unit {
                         coef: One::one(),
                         xpow: q.power(&BigInt::from(2)),
                         ypow: BigInt::from(0),
                       }.to_pol() * division_polynomial::psi(a, b, &ql).power(&BigInt::from(2)));
        //println!("n2:{} {}", line!(), n2.to_string());

        let n3 = (division_polynomial::phi(a, b, &ql)
                     - Unit {
                         coef: One::one(),
                         xpow: q.power(&BigInt::from(2)),
                         ypow: BigInt::from(0),
                       }.to_pol() * division_polynomial::psi(a, b, &ql).power(&BigInt::from(2))).power(&BigInt::from(2));
        //println!("n3:{} {}", line!(), n3.to_string());

        let b_num1 = n1 + n2 * n3; 
        println!("b_num1:{} {}", line!(), b_num1.to_string());

        let b_num1 = b_num1.modular(q);
        println!("b_num1:{} {}", line!(), b_num1.to_string());

        let b_num1 = b_num1.ec_reduction(a, b);

        let b_num1 = b_num1.modular(q);

        let b_den1 = (division_polynomial::psi(a, b, &ql).power(&BigInt::from(2)))
                 * ((division_polynomial::phi(a, b, &ql) - 
                        Unit {
                            coef: BigInt::from(1),
                            xpow: q.power(&BigInt::from(2)),
                            ypow: BigInt::from(0),
                        }.to_pol()
                   * division_polynomial::psi(a, b, &ql).power(&BigInt::from(2))).power(&BigInt::from(2)));
        let b_den1 = b_den1.modular(q);
        
        println!("b_den1:{} {}", line!(), b_den1.to_string());

        let mut found = false;
        let mut jj: BigInt = BigInt::from(0);
        for j in num_iter::range(BigInt::from(1), jmax + BigInt::from(1)) {
            println!("{} j:{}", line!(), j.to_string());
            // x
            let num2 = division_polynomial::phi(a, b, &j).modular(q).to_frob(q);
            println!("num2:{} {}", line!(), num2.to_string());

            let den2 = division_polynomial::psi(a, b, &j).modular(q).to_frob(q).power(&BigInt::from(2));
            println!("den2:{} {}", line!(), den2.to_string());

            let p1 = &b_num1 * &den2 - &num2 * &b_den1;
            let p1 = p1.modular(q);

            let p1 = p1.ec_reduction(&a, b);
            //println!("p1:{} {}", line!(), p1.to_string());

            let p1 = p1.modular(q);
            let psil = division_polynomial::psi(a, b, &l).modular(q);
            println!("{} psi({}):{}", line!(), l, psil.to_string());

            let p1 = p1.polynomial_modular(&psil, q);

            let p1 = &p1.modular(q);
            println!("p1:{} {}", line!(), p1.to_string());

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
            let d = (division_polynomial::omega(a, b, &ql) - 
                    Unit {
                        coef: BigInt::from(1),
                        xpow: BigInt::from(0),
                        ypow: q.power(&BigInt::from(2)),
                    }.to_pol()
                    ) * division_polynomial::psi(a, b, &ql).power(&BigInt::from(3));
            let d = d.modular(q);
            //println!("d:{} {}", line!(), d.to_string());
            let e = - division_polynomial::phi(a, b, &ql) +
                    Unit {
                        coef: BigInt::from(2),
                        xpow: q.power(&BigInt::from(2)),
                        ypow: BigInt::from(0),
                    }.to_pol() * division_polynomial::psi(a, b, &ql).power(&BigInt::from(2));
            let e = e.modular(q);
            //println!("e:{} {}", line!(), e.to_string());

            let e = e.ec_reduction(a, b);
            //println!("e:{} {}", line!(), e.to_string());

            let f = division_polynomial::phi(a, b, &ql) - Unit {
                    coef: BigInt::from(1),
                    xpow: q.power(&BigInt::from(2)),
                    ypow: BigInt::from(0),
                }.to_pol() * division_polynomial::psi(a, b, &ql).power(&BigInt::from(2));
            let f = f.modular(q);
            //println!("f:{} {}", line!(), f.to_string());

            let f = f.ec_reduction(a, b);
            //println!("f:{} {}", line!(), f.to_string());

            let num1 = &d * (e * f.power(&BigInt::from(2)) - d.clone().power(&BigInt::from(2))) -
                Unit {
                    coef: BigInt::from(1),
                    xpow: BigInt::from(0),
                    ypow: q.power(&BigInt::from(2)),
                }.to_pol() * 
                division_polynomial::psi(a, b, &ql).power(&BigInt::from(3)) * f.power(&BigInt::from(3));
            let num1 = num1.modular(q);
            println!("num1:{} {}", line!(), num1.to_string());

            let num1 = num1 * division_polynomial::psi(a, b, &ql);
            let num1 = num1.modular(q);
            println!("num1:{} {}", line!(), num1.to_string());

            let den1 = division_polynomial::psi(a, b, &ql).power(&BigInt::from(3)) * f.power(&BigInt::from(3)) * division_polynomial::psi(a, b, &ql);
            let den1 = den1.modular(q);
            println!("den1:{} {}", line!(), den1.to_string());

            let num2 = division_polynomial::omega(a, b, &jj).to_frob(q);
            let num2 = num2.modular(q);
            println!("num2:{} {}", line!(), num2.to_string());

            let den2 = division_polynomial::psi(a, b, &jj).to_frob(q).power(&BigInt::from(3));
            println!("den2:{} {}", line!(), den2.to_string());

            let p1 = num1 * den2 - num2 * den1;
            let p2 = p1.ec_reduction(a, b);
            println!("p2:{} {}", line!(), p2.to_string());

            let p3: Polynomial;
            if p2.has_y() {
                p3 = p2 / UnitBuilder::new().coef(1).ypow(1).finalize().to_pol();
            } else {
                p3 = p2.clone();
            }
            println!("{} p3:{}", line!(), p3.to_string());

            let p4 = p3.modular(q);
            let psil = division_polynomial::psi(a, b, &l).modular(q);
            println!("psi({}):{}", l.to_string(), psil.to_string());

            let p6 = p4.polynomial_modular(&psil, q);
            println!("pol % psi({}):{}", l.to_string(), p6.to_string());

            if p6.is_zero() {
                println!("a = {} mod {}", &jj.to_string(), l.to_string());
            } else {
                println!("a = {} mod {}", (-jj.clone()).to_string(), l.to_string());
            }
        } else {
            // # (d)
            println!("x not found");
            let mut found = false;
            let mut w = BigInt::from(1);
            for i in num_iter::range(BigInt::from(1), l.clone()) {
                if i.power(&BigInt::from(2)).rem_floor(l) == ql {
                    found = true;
                    w = i.clone();
                    break;
                }
            }
            if !found {
                println!("a = 0 mod l because of w not found (d)");
            } else {
                // (e) x
                let p1 = Unit {
                        coef: BigInt::from(1),
                        xpow: q.clone(),
                        ypow: BigInt::from(0),
                    }.to_pol() * division_polynomial::psi(a, b, &w).power(&BigInt::from(2)) - division_polynomial::phi(a, b, &w);
                let p1 = p1.modular(q);
                let p2 = p1.ec_reduction(a, b);
                println!("{} p2:{}", line!(), p2.to_string());
                let p3 = division_polynomial::psi(q, b, &l).modular(q);
                println!("{} p3:{}", line!(), p3.to_string());

                let p4 = p2.polynomial_modular(&p3, q);
                println!("{} p4 {}", line!(), p4.to_string());
                if !p4.is_zero() {
                    println!("a = 0 mod {} because of gcd = 1 (e)", l.to_string());
                } else {
                    // (e) y
                    let p5 = (Unit {
                                    coef: BigInt::from(1),
                                    xpow: BigInt::from(0),
                                    ypow: q.clone(),
                                }.to_pol() * division_polynomial::psi(a, b, &w).power(&BigInt::from(3)) - division_polynomial::omega(a, b, &w)) /
                            UnitBuilder::new().coef(1).ypow(1).finalize().to_pol();
                    let p5 = p5.modular(q);
                    let p6 = p5.ec_reduction(a, b);
                    println!("{}", p6.to_string());
                    if p6.is_gcd_one(&division_polynomial::psi(a, b, &l), q) {
                        println!("a = {} mod {}", (- BigInt::from(2) * w).to_string(), l.to_string());
                    } else {
                        println!("a = {} mod {}", (BigInt::from(2) * w).to_string(), l.to_string());
                    }
                }
            }
        }
    }
}

#[test]
fn schoof_test() {
    schoof(&BigInt::from(2), &BigInt::from(1), &BigInt::from(19));
}

