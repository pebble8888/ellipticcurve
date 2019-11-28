use primes;
use num_integer::Integer;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;
use super::bigint::Power;
use super::j_invariant;
use super::subscripted_variable;
use std::ops::Neg;

/// calculate modular polynomial
pub fn subscripted_variable_modular_polynomial(p: u64) -> polynomial::Polynomial {
    assert!(p >= 2);
    assert!(primes::is_prime(p));
    let mut pol = term_builder::TermBuilder::new()
                .xpow(&BigInt::from(p+1))
                .build().to_pol();
    pol += term_builder::TermBuilder::new()
                .ypow(&BigInt::from(p+1))
                .build().to_pol();

    for i in num_iter::range(0, p+1) {
        pol += term_builder::TermBuilder::new()
            .xpow(&BigInt::from(i))
            .ypow(&BigInt::from(i))
            .variable_ij(i, i)
            .build().to_pol();
    }
    for i in num_iter::range(0, p+1) {
        for j in num_iter::range(i+1, p+1) {
            pol += term_builder::TermBuilder::new()
                .xpow(&BigInt::from(i))
                .ypow(&BigInt::from(j))
                .variable_ij(i, j)
                .build().to_pol();
            pol += term_builder::TermBuilder::new()
                .xpow(&BigInt::from(j))
                .ypow(&BigInt::from(i))
                .variable_ij(i, j)
                .build().to_pol();
        }
    }
    pol
}

pub fn subscripted_variable_modular_polynomial_list(p: u64) -> Vec<polynomial::Polynomial> {
    let pol_q = subscripted_variable_modular_polynomial_q(p);
    let pp = BigInt::from(p);
    let min = BigInt::from(- pp.power(2) - pp);
    let mut v: Vec<polynomial::Polynomial> = Vec::new();
    for i in num_iter::range(min, One::one()) {
        v.push(pol_q.to_q_power_coef(&i));
    }
    v
}

pub fn subscripted_variable_modular_polynomial_q(p: u64) -> polynomial::Polynomial {
    let pol = subscripted_variable_modular_polynomial(p);
    // j-invariant q order needs p^2 + p
    let j_order = p * p + p;
    let j = j_invariant::j_invariant(j_order);
    let j_q_power = j.to_q_power(p as i64);
    let pol1 = pol.eval_x_polynomial(&j);
    let pol2 = pol1.eval_y_polynomial(&j_q_power);
    let pol_q = pol2.omit_high_order_q(0); 
    pol_q
}

pub fn modular_polynomial(p: u64) -> polynomial::Polynomial {
    // TODO: solver simultanous equation1 to another file
    let list = subscripted_variable_modular_polynomial_list(p);
    // (i) LU
    let converter = subscripted_variable::SubscriptedVariableConverter::new(p);
    let row_count: usize = list.len();
    let col_count: usize = converter.count() as usize + 1;
    let mut a: Vec<Vec<BigInt>> = vec![vec![BigInt::from(0); col_count]; row_count];
    for row in num_iter::range(0, row_count) {
        for col in num_iter::range(0, col_count) {
            if col < col_count - 1 {
                let variable = converter.variable_from_index(col as u64);
                let coef = list[row].to_variable_coef(variable);
                a[row][col] = coef;
            } else if col == col_count - 1 {
                let coef = list[row].to_variable_coef(subscripted_variable::SubscriptedVariable::new());
                a[row][col] = coef;
            } else {
                assert!(false);
            }
        }
    }

    let mut row: Vec<usize> = Vec::new();
    for i in num_iter::range(0, row_count) {
        row.push(i);
    }

    for i in num_iter::range(0, col_count - 1) {
        // (I)
        for j in num_iter::range(i, row_count) {
            if a[row[j]][i] != BigInt::from(0) {
                if i != j {
                    // row swap
                    let b = row[i];
                    row[i] = row[j];
                    row[j] = b;
                }
                break;
            }
        }
        let c1 = a[row[i]][i].clone();
        assert!(c1 != BigInt::from(0));
        // (II)
        for j in num_iter::range(i+1, row_count) {
            let c2 = a[row[j]][i].clone();
            if c2 != BigInt::from(0) {
                let lcm = Integer::lcm(&c1, &c2); 
                let cc1 = lcm.clone() / c1.clone();
                //assert_eq!(cc1.clone() * c1.clone(), lcm.clone());
                let cc2 = lcm.clone() / c2.clone();
                //assert_eq!(cc2.clone() * c2.clone(), lcm.clone());
                for k in num_iter::range(0, col_count) {
                    a[row[j]][k] *= cc2.clone();
                }
                for k in num_iter::range(0, col_count) {
                    let val = a[row[i]][k].clone() * cc1.clone();
                    a[row[j]][k] -= val;
                }
            }
        }
    }

    // (ii) Substitute
    let col_prim = col_count -1;
    for i in num_iter::range(1, col_count - 1).rev() {
        // diagonal coef
        let diag = a[row[i]][i].clone();
        let prim = a[row[i]][col_prim].clone();
        //assert_eq!(prim.mod_floor(&diag.clone()), BigInt::from(0));
        if diag != One::one() {
            a[row[i]][i] = 1.into();
            a[row[i]][col_prim] = prim / diag;
        }
        let val = a[row[i]][col_prim].clone();
        for j in num_iter::range(0, i) {
            let t_val = a[row[j]][i].clone();
            a[row[j]][col_prim] -= t_val * val.clone();
            a[row[j]][i] = Zero::zero();
        }
    }

    let mut pol = polynomial::Polynomial::new();
    for i in num_iter::range(0, col_count - 1) {
        let variable = converter.variable_from_index(i as u64);
        let val = a[row[i]][col_prim].clone();
        pol += term_builder::TermBuilder::new()
            .coef(&val.clone().neg())
            .xpow(variable.i as i64)
            .ypow(variable.j as i64)
            .build();
        if variable.i != variable.j {
            pol += term_builder::TermBuilder::new()
                .coef(&val.clone().neg())
                .xpow(variable.j as i64)
                .ypow(variable.i as i64)
                .build();
        }
    }
    pol += term_builder::TermBuilder::new().xpow((p as i64)+1).build();
    pol += term_builder::TermBuilder::new().ypow((p as i64)+1).build();
    return pol;
}

#[test]
fn modular_polynomial_test2() {
    let pol = modular_polynomial(2);
    assert_eq_str!(pol, "x^3 - x^2 y^2 + 1488 x^2 y - 162000 x^2 + 1488 x y^2 + 40773375 x y + 8748000000 x + y^3 - 162000 y^2 + 8748000000 y - 157464000000000");
}

#[test]
fn modular_polynomial_test3() {
    let pol = modular_polynomial(3);
    assert_eq_str!(pol, "x^4 - x^3 y^3 + 2232 x^3 y^2 - 1069956 x^3 y + 36864000 x^3 + 2232 x^2 y^3 + 2587918086 x^2 y^2 + 8900222976000 x^2 y + 452984832000000 x^2 - 1069956 x y^3 + 8900222976000 x y^2 - 770845966336000000 x y + 1855425871872000000000 x + y^4 + 36864000 y^3 + 452984832000000 y^2 + 1855425871872000000000 y");
}

#[test]
fn modular_polynomial_test5() {
    let pol = modular_polynomial(5);
    assert_eq_str!(pol, "x^6 - x^5 y^5 + 3720 x^5 y^4 - 4550940 x^5 y^3 + 2028551200 x^5 y^2 - 246683410950 x^5 y + 1963211489280 x^5 + 3720 x^4 y^5 + 1665999364600 x^4 y^4 + 107878928185336800 x^4 y^3 + 383083609779811215375 x^4 y^2 + 128541798906828816384000 x^4 y + 1284733132841424456253440 x^4 - 4550940 x^3 y^5 + 107878928185336800 x^3 y^4 - 441206965512914835246100 x^3 y^3 + 26898488858380731577417728000 x^3 y^2 - 192457934618928299655108231168000 x^3 y + 280244777828439527804321565297868800 x^3 + 2028551200 x^2 y^5 + 383083609779811215375 x^2 y^4 + 26898488858380731577417728000 x^2 y^3 + 5110941777552418083110765199360000 x^2 y^2 + 36554736583949629295706472332656640000 x^2 y + 6692500042627997708487149415015068467200 x^2 - 246683410950 x y^5 + 128541798906828816384000 x y^4 - 192457934618928299655108231168000 x y^3 + 36554736583949629295706472332656640000 x y^2 - 264073457076620596259715790247978782949376 x y + 53274330803424425450420160273356509151232000 x + y^6 + 1963211489280 y^5 + 1284733132841424456253440 y^4 + 280244777828439527804321565297868800 y^3 + 6692500042627997708487149415015068467200 y^2 + 53274330803424425450420160273356509151232000 y + 141359947154721358697753474691071362751004672000");
}

#[test]
fn subscripted_variable_modular_polynomial_p2_test() {
    let p = 2;
    let pol = subscripted_variable_modular_polynomial(p);
    assert_eq_str!(pol, "x^3 + x^2 y^2 c_2_2 + x^2 y c_1_2 + x^2 c_0_2 + x y^2 c_1_2 + x y c_1_1 + x c_0_1 + y^3 + y^2 c_0_2 + y c_0_1 + c_0_0");

    let pol_q = subscripted_variable_modular_polynomial_q(p);
    assert_eq_str!(pol_q, "941847590656704 c_2_2 + 126112980648 c_1_2 + 22047296 c_1_1 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 + 2710404480 + 10291429500960 q^-1 c_2_2 + 1495268650 q^-1 c_1_2 + 197628 q^-1 c_1_1 + 1488 q^-1 c_0_2 + q^-1 c_0_1 + 2251260 q^-1 + 73885159932 q^-2 c_2_2 + 23548880 q^-2 c_1_2 + 744 q^-2 c_1_1 + 1489 q^-2 c_0_2 + q^-2 c_0_1 + 2253492 q^-2 + 338165056 q^-3 c_2_2 + 199860 q^-3 c_1_2 + q^-3 c_1_1 + q^-3 + 948792 q^-4 c_2_2 + 745 q^-4 c_1_2 + q^-4 c_0_2 + 2232 q^-4 + 1488 q^-5 c_2_2 + q^-5 c_1_2 + q^-6 c_2_2 + q^-6");

    let list = subscripted_variable_modular_polynomial_list(p);
    assert_eq!(list.len(), 7);
    assert_eq_str!(list[0], "c_2_2 + 1");
    assert_eq_str!(list[1], "1488 c_2_2 + c_1_2");
    assert_eq_str!(list[2], "948792 c_2_2 + 745 c_1_2 + c_0_2 + 2232");
    assert_eq_str!(list[3], "338165056 c_2_2 + 199860 c_1_2 + c_1_1 + 1");
    assert_eq_str!(list[4], "73885159932 c_2_2 + 23548880 c_1_2 + 744 c_1_1 + 1489 c_0_2 + c_0_1 + 2253492");
    assert_eq_str!(list[5], "10291429500960 c_2_2 + 1495268650 c_1_2 + 197628 c_1_1 + 1488 c_0_2 + c_0_1 + 2251260");
    assert_eq_str!(list[6], "941847590656704 c_2_2 + 126112980648 c_1_2 + 22047296 c_1_1 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 + 2710404480");
}

#[test]
fn subscripted_variable_modular_polynomial_p3_test() {
    let p = 3;
    let pol = subscripted_variable_modular_polynomial(p);
    assert_eq_str!(pol, "x^4 + x^3 y^3 c_3_3 + x^3 y^2 c_2_3 + x^3 y c_1_3 + x^3 c_0_3 + x^2 y^3 c_2_3 + x^2 y^2 c_2_2 + x^2 y c_1_2 + x^2 c_0_2 + x y^3 c_1_3 + x y^2 c_1_2 + x y c_1_1 + x c_0_1 + y^4 + y^3 c_0_3 + y^2 c_0_2 + y c_0_1 + c_0_0");

    let pol_q = subscripted_variable_modular_polynomial_q(p);
    assert_eq_str!(pol_q, "4197391688509451809042817340 c_3_3 + 65201219062152355627872 c_2_3 + 1769523054225123330 c_2_2 + 45070853378831784 c_1_3 + 15329636199360 c_1_2 + 864853506 c_1_1 + 2710404480 c_0_3 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 + 4084248062160 + 109642480943384830620762360 q^-1 c_3_3 + 1768996792073642463801 q^-1 c_2_3 + 45738435798911040 q^-1 c_2_2 + 1345112306562240 q^-1 c_1_3 + 437662034132 q^-1 c_1_2 + 21494504 q^-1 c_1_1 + 2251260 q^-1 c_0_3 + 1488 q^-1 c_0_2 + q^-1 c_0_1 + 3491078528 q^-1 + 2319314637883533604808364 q^-2 c_3_3 + 57383938307316236432 q^-2 c_2_3 + 833607524819048 q^-2 c_2_2 + 90830762088165 q^-2 c_1_3 + 20874771304 q^-2 c_1_2 + 196884 q^-2 c_1_1 + 2232 q^-2 c_0_3 + q^-2 c_0_2 + 4108752 q^-2 + 38360228325195786383216 q^-3 c_3_3 + 1807266275628309954 q^-3 c_2_3 + 9791534543904 q^-3 c_2_2 + 6184170973560 q^-3 c_1_3 + 866354346 q^-3 c_1_2 + 744 q^-3 c_1_1 + 2251261 q^-3 c_0_3 + 1488 q^-3 c_0_2 + q^-3 c_0_1 + 3491081504 q^-3 + 480257921200323487545 q^-4 c_3_3 + 45943884263343552 q^-4 c_2_3 + 72476838420 q^-4 c_2_2 + 381181215440 q^-4 c_1_3 + 21496736 q^-4 c_1_2 + q^-4 c_1_1 + q^-4 + 4461984465203723760 q^-5 c_3_3 + 834399255041138 q^-5 c_2_3 + 335952400 q^-5 c_2_2 + 20685303576 q^-5 c_1_3 + 196885 q^-5 c_1_2 + 30459141465291828 q^-6 c_3_3 + 9793594541808 q^-6 c_2_3 + 947304 q^-6 c_2_2 + 865960579 q^-6 c_1_3 + 744 q^-6 c_1_2 + 2232 q^-6 c_0_3 + q^-6 c_0_2 + 4108752 q^-6 + 151527078622080 q^-7 c_3_3 + 72480196752 q^-7 c_2_3 + 1488 q^-7 c_2_2 + 21495992 q^-7 c_1_3 + q^-7 c_1_2 + 541783100214 q^-8 c_3_3 + 335955376 q^-8 c_2_3 + q^-8 c_2_2 + 196884 q^-8 c_1_3 + 1355204472 q^-9 c_3_3 + 947305 q^-9 c_2_3 + 744 q^-9 c_1_3 + q^-9 c_0_3 + 2976 q^-9 + 2251260 q^-10 c_3_3 + 1488 q^-10 c_2_3 + q^-10 c_1_3 + 2232 q^-11 c_3_3 + q^-11 c_2_3 + q^-12 c_3_3 + q^-12");

    let list = subscripted_variable_modular_polynomial_list(p);
    assert_eq!(list.len(), 13);
    assert_eq_str!(list[0], "c_3_3 + 1");
    assert_eq_str!(list[1], "2232 c_3_3 + c_2_3");
    assert_eq_str!(list[2], "2251260 c_3_3 + 1488 c_2_3 + c_1_3");
    assert_eq_str!(list[3], "1355204472 c_3_3 + 947305 c_2_3 + 744 c_1_3 + c_0_3 + 2976");
    assert_eq_str!(list[4], "541783100214 c_3_3 + 335955376 c_2_3 + c_2_2 + 196884 c_1_3");
    assert_eq_str!(list[5], "151527078622080 c_3_3 + 72480196752 c_2_3 + 1488 c_2_2 + 21495992 c_1_3 + c_1_2");
    assert_eq_str!(list[6], "30459141465291828 c_3_3 + 9793594541808 c_2_3 + 947304 c_2_2 + 865960579 c_1_3 + 744 c_1_2 + 2232 c_0_3 + c_0_2 + 4108752");
    assert_eq_str!(list[7], "4461984465203723760 c_3_3 + 834399255041138 c_2_3 + 335952400 c_2_2 + 20685303576 c_1_3 + 196885 c_1_2");
    assert_eq_str!(list[8], "480257921200323487545 c_3_3 + 45943884263343552 c_2_3 + 72476838420 c_2_2 + 381181215440 c_1_3 + 21496736 c_1_2 + c_1_1 + 1");
    assert_eq_str!(list[9], "38360228325195786383216 c_3_3 + 1807266275628309954 c_2_3 + 9791534543904 c_2_2 + 6184170973560 c_1_3 + 866354346 c_1_2 + 744 c_1_1 + 2251261 c_0_3 + 1488 c_0_2 + c_0_1 + 3491081504");
    assert_eq_str!(list[10], "2319314637883533604808364 c_3_3 + 57383938307316236432 c_2_3 + 833607524819048 c_2_2 + 90830762088165 c_1_3 + 20874771304 c_1_2 + 196884 c_1_1 + 2232 c_0_3 + c_0_2 + 4108752");
    assert_eq_str!(list[11], "109642480943384830620762360 c_3_3 + 1768996792073642463801 c_2_3 + 45738435798911040 c_2_2 + 1345112306562240 c_1_3 + 437662034132 c_1_2 + 21494504 c_1_1 + 2251260 c_0_3 + 1488 c_0_2 + c_0_1 + 3491078528");
    assert_eq_str!(list[12], "4197391688509451809042817340 c_3_3 + 65201219062152355627872 c_2_3 + 1769523054225123330 c_2_2 + 45070853378831784 c_1_3 + 15329636199360 c_1_2 + 864853506 c_1_1 + 2710404480 c_0_3 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 + 4084248062160");
}

/*
#[test]
fn modular_polynomial_test2() {
    let pol = modular_polynomial(2);
    assert_eq_str!(pol, "x^3 - x^2 y^2 + 1488 x^2 y - 162000 x^2 + 1488 x y^2 + 40773375 x y + 8748000000 x + y^3 - 162000 y^2 + 8748000000 y - 157464000000000");
}

#[test]
fn modular_polynoial_test3() {
    let pol = modular_polynomial(3);
    assert_eq_str!(pol, "x^4 - x^3 y^3 + 2232 x^3 y^2 - 1069956 x^3 y + 36864000 x^3 + 2232 x^2 y^3 + 2587918086 x^2 y^2 + 8900222976000 x^2 y + 452984832000000 x^2 - 1069956 x y^3 + 8900222976000 x y^2 - 770845966336000000 x y + 1855425871872000000000 x + y^4 + 36864000 y^3 + 452984832000000 y^2 + 1855425871872000000000 y");
}
*/

