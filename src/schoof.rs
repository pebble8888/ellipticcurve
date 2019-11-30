extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;
use num_traits::{Zero, One};
use num_integer::Integer;
use num_bigint::BigInt;
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;
use super::division_polynomial;
use crate::bigint;
use crate::bigint::{Power};

type TermBuilder = term_builder::TermBuilder;
type Polynomial = polynomial::Polynomial;

/// Calculate the order of the Elliptic curve y^2 = x^3 + a x + b
/// schoof algorithm
pub fn schoof(a: &BigInt, b: &BigInt, q: &BigInt) -> Vec<bigint::ModResult> {
    let mut mod_result: Vec<bigint::ModResult> = Vec::new();

    // l = 2
    let l: BigInt = 2.into();
    println!("{} l:{}", line!(), l);
    let pol_l2 = TermBuilder::new().xpow(q).build() - TermBuilder::new().coef(1).xpow(1).build();
    let pol_standard = TermBuilder::new().xpow(3).build() + 
                TermBuilder::new().coef(a).xpow(1).build() +
                TermBuilder::new().coef(b).build();
    let j2: BigInt;
    if pol_standard.is_gcd_one(&pol_l2, q) {
        // no common root
        j2 = One::one();
    } else {
        // has common root
        j2 = Zero::zero();
    }
    mod_result.push(bigint::ModResult { l: l.clone(), r: j2.clone() });
    println!("a = {} mod 2", j2); 

    // l >= 3
    // TODO:
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

        let nn2 = TermBuilder::new().ypow(&q.power(2)).build()
                * division_polynomial::psi(a, b, &ql).power(3);

        let nn2 = nn2.reduction_modular(a, b, q);

        let mut n1 = (nn1 - nn2).power(2);

        n1.modular_assign(q);

        let n2 = - (division_polynomial::phi(a, b, &ql)
                     + TermBuilder::new().xpow(&q.power(2)).build()
                       * division_polynomial::psi(a, b, &ql).power(2));


        let n3 = (division_polynomial::phi(a, b, &ql)
                     - TermBuilder::new().xpow(&q.power(2)).build()
                       * division_polynomial::psi(a, b, &ql).power(2)).power(2);

        let mut num1 = n1 + n2 * n3;
        num1.modular_assign(q);
        let num1 = num1.reduction_modular(a, b, q);

        let den1 = (division_polynomial::psi(a, b, &ql).power(2))
                 * ((division_polynomial::phi(a, b, &ql) - 
                        TermBuilder::new().xpow(&q.power(2)).build()
                   * division_polynomial::psi(a, b, &ql).power(2)).power(2));
        
        let den1 = den1.reduction_modular(a, b, q);

        let psi_l = division_polynomial::psi(a, b, &l).modulo(q);
        println!("{} psi({}):{}", line!(), l, psi_l);

        let mut found = false;
        let mut jj: BigInt = Zero::zero();
        for j in num_iter::range(BigInt::from(1), jmax + BigInt::from(1)) {
            println!("{} j:{}", line!(), j);
            // x
            let num2 = division_polynomial::phi(a, b, &j).to_frob(q);
            let num2 = num2.reduction_modular(a, b, q);

            let den2 = division_polynomial::psi(a, b, &j).to_frob(q).power(2);

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
            let g = TermBuilder::new().ypow(&q.power(2)).build().to_pol();
            let g = g.reduction_modular(a, b, q);

            let mut omg = division_polynomial::omega(a, b, &ql);
            omg.modular_assign(q);

            let d = omg - &g * division_polynomial::psi(a, b, &ql).power(3);
            let d = d.reduction_modular(a, b, q);

            let mut e = - division_polynomial::phi(a, b, &ql) +
                    TermBuilder::new().coef(2).xpow(&q.power(2)).build()
                    * division_polynomial::psi(a, b, &ql).power(2);
            e.modular_assign(q);

            let e = e.reduction_modular(a, b, q);

            let mut f = division_polynomial::phi(a, b, &ql) - 
                TermBuilder::new().xpow(&q.power(2)).build() * division_polynomial::psi(a, b, &ql).power(2);
            f.modular_assign(q);

            let f = f.reduction_modular(a, b, q);

            let num3 = &d * (e * f.power(2) - d.clone().power(2)) - &g * division_polynomial::psi(a, b, &ql).power(3) * f.power(3);
            let mut num3 = num3.reduction_modular(a, b, q);

            num3 *= division_polynomial::psi(a, b, &ql);
            let num3 = num3.reduction_modular(a, b, q);

            let den3 = division_polynomial::psi(a, b, &ql).power(3) * f.power(3) * division_polynomial::psi(a, b, &ql);
            let den3 = den3.reduction_modular(a, b, q);

            let num4 = division_polynomial::omega(a, b, &jj).to_frob(q);
            let num4 = num4.reduction_modular(a, b, q);

            let den4 = division_polynomial::psi(a, b, &jj).to_frob(q).power(3);
            let den4 = den4.reduction_modular(a, b, q);

            let mut p7 = num3 * den4 - num4 * den3;
            p7.modular_assign(q);

            let p8 = p7.reduction_modular(a, b, q);

            let mut p9: Polynomial;
            if p8.has_y() {
                p9 = p8 / TermBuilder::new().ypow(1).build();
            } else {
                p9 = p8.clone();
            }

            p9.modular_assign(q);

            let p9 = p9.polynomial_modular(&psi_l, q);
            println!("{} pol % psi({}):{}", line!(), l, p9);

            if !p9.is_zero() {
                jj = -jj;
            }
            if jj < BigInt::from(0) {
                jj += l;
            }
            mod_result.push(bigint::ModResult { l: l.clone(), r: jj.clone() });
            println!("(iii) y  a = {} mod {}", jj.clone(), l);
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
                mod_result.push(bigint::ModResult { l: l.clone(), r: BigInt::from(0) });
                println!("(d)  a = 0 mod l because of w not found (d)");
            } else {
                // (e) x
                let mut p12 = TermBuilder::new().xpow(&q.clone()).build()
                       * division_polynomial::psi(a, b, &w).power(2) - division_polynomial::phi(a, b, &w);
                p12.modular_assign(q);
                let p13 = p12.reduction_modular(a, b, q);

                let p14 = division_polynomial::psi(q, b, &l).modulo(q);

                let p15 = p13.polynomial_modular(&p14, q);

                if !p15.is_zero() {
                    mod_result.push(bigint::ModResult { l: l.clone(), r: BigInt::from(0) });
                    println!("(e) y  a = 0 mod {}", l);
                } else {
                    // (e) y
                    let mut p16 = (TermBuilder::new().ypow(&a.clone()).build()
                              * division_polynomial::psi(a, b, &w).power(3) - division_polynomial::omega(a, b, &w)) /
                            TermBuilder::new().ypow(1).build();
                    p16.modular_assign(q);
                    let p17 = p16.reduction_modular(a, b, q);
                    let mut ww: BigInt;
                    if p17.is_gcd_one(&division_polynomial::psi(a, b, &l), q) {
                        ww = - BigInt::from(2) * w;
                    } else {
                        ww = BigInt::from(2) * w;
                    }
                    if ww < BigInt::from(0) {
                        ww += l;
                    }
                    mod_result.push(bigint::ModResult { l: l.clone(), r: ww.clone() });
                    println!("(e) y  a = {} mod {}", &ww, l);
                }
            }
        }
    }
    mod_result
}

#[test]
#[ignore]
fn schoof_test7() {
    use crate::bigint::chinese_remainder;

    let q = 7;
    let mod_result = schoof(&BigInt::from(2), &BigInt::from(1), &BigInt::from(q));
    assert_eq!(mod_result.len(), 3);
    assert_eq!(mod_result[0].l, BigInt::from(2));
    assert_eq!(mod_result[0].r, BigInt::from(1));
    assert_eq!(mod_result[1].l, BigInt::from(3));
    assert_eq!(mod_result[1].r, BigInt::from(0));
    assert_eq!(mod_result[2].l, BigInt::from(5));
    assert_eq!(mod_result[2].r, BigInt::from(3));

    let result = chinese_remainder(&mod_result);
    assert_eq!(result.r, BigInt::from(3));
    assert_eq!(result.l, BigInt::from(30)); // 2 * 3 * 5
}

#[test]
#[ignore]
fn schoof_test19() {
    use crate::bigint::chinese_remainder;
    let q = 19;
    let mod_result = schoof(&BigInt::from(2), &BigInt::from(1), &BigInt::from(q));
    assert_eq!(mod_result.len(), 3);
    assert_eq!(mod_result[0].l, BigInt::from(2));
    assert_eq!(mod_result[0].r, BigInt::from(1));
    assert_eq!(mod_result[1].l, BigInt::from(3));
    assert_eq!(mod_result[1].r, BigInt::from(2));
    assert_eq!(mod_result[2].l, BigInt::from(5));
    assert_eq!(mod_result[2].r, BigInt::from(3));

    let result = chinese_remainder(&mod_result);
    assert_eq!(result.r, BigInt::from(23));
    assert_eq!(result.l, BigInt::from(30)); // 2 * 3 * 5
}

