use primes;
use num_bigint::BigInt;
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;
use super::bigint::Power;
use super::j_invariant;

type TermBuilder = term_builder::TermBuilder;
type Polynomial = polynomial::Polynomial;

/// create modular polynomial
pub fn modular_polynomial(n: u64) -> Polynomial {
    assert!(n >= 2);
    if n == 2 {
          TermBuilder::new().coef(1).xpow(3).build()
        + TermBuilder::new().coef(1).ypow(3).build()
        + TermBuilder::new().coef(-162_000).xpow(2).build()
        + TermBuilder::new().coef(-162_000).ypow(2).build()
        + TermBuilder::new().coef(1_488).xpow(2).ypow(1).build()
        + TermBuilder::new().coef(1_488).ypow(2).xpow(1).build()
        + TermBuilder::new().coef(-1).xpow(2).ypow(2).build()
        + TermBuilder::new().coef(8_748_000_000).xpow(1).build()
        + TermBuilder::new().coef(8_748_000_000).ypow(1).build()
        + TermBuilder::new().coef(40_773_375).xpow(1).ypow(1).build()
        + TermBuilder::new().coef(-157_464_000_000_000).build()
    } else if n == 3 {
        let c10 = BigInt::from(1_855_425_871_872_000_000_000_i128);
        let pol =
          TermBuilder::new().coef(1).xpow(4).build()
        + TermBuilder::new().coef(1).ypow(4).build()
        + TermBuilder::new().coef(-1).xpow(3).ypow(3).build()
        + TermBuilder::new().coef(2_232).xpow(3).ypow(2).build()
        + TermBuilder::new().coef(2_232).ypow(3).xpow(2).build()
        + TermBuilder::new().coef(-1_069_956).xpow(3).ypow(1).build()
        + TermBuilder::new().coef(-1_069_956).xpow(1).ypow(3).build()
        + TermBuilder::new().coef(36_864_000).xpow(3).build()
        + TermBuilder::new().coef(36_864_000).ypow(3).build()
        + TermBuilder::new().coef(2_587_918_086).xpow(2).ypow(2).build()
        + TermBuilder::new().coef(8_900_222_976_000).xpow(2).ypow(1).build()
        + TermBuilder::new().coef(8_900_222_976_000).ypow(2).xpow(1).build()
        + TermBuilder::new().coef(452_984_832_000_000).xpow(2).build()
        + TermBuilder::new().coef(452_984_832_000_000).ypow(2).build()
        + TermBuilder::new().coef(-770_845_966_336_000_000).xpow(1).ypow(1).build()
        + TermBuilder::new().coef(&c10).xpow(1).build()
        + TermBuilder::new().coef(&c10).ypow(1).build();
        return pol;
    } else if n == 5 {
        let b = BigInt::from(1_000_000_000_000_000_000_000_000_000_000_000_000_i128);
        let c00 = BigInt::from(141_359_947_154_i64) * b.power(1)
                + BigInt::from(721_358_697_753_474_691_071_362_751_004_672_000_i128);
        let c10 = BigInt::from(53_274_330) * b.power(1)
                + BigInt::from(803_424_425_450_420_160_273_356_509_151_232_000_i128);
        let c11 = - BigInt::from(264_073) * b.power(1)
                - BigInt::from(457_076_620_596_259_715_790_247_978_782_949_376_i128);
        let c20 = BigInt::from(6_692) * b.power(1)
                + BigInt::from(500_042_627_997_708_487_149_415_015_068_467_200_i128);
        let c21 = BigInt::from(36) * b.power(1)
                + BigInt::from(554_736_583_949_629_295_706_472_332_656_640_000_i128);
        let c22 = BigInt::from(5_110_941_777_552_418_083_110_765_199_360_000_i128);
        let c30 = BigInt::from(280_244_777_828_439_527_804_321_565_297_868_800_i128);
        let c31 = - BigInt::from(192_457_934_618_928_299_655_108_231_168_000_i128);
        let c32 = BigInt::from(26_898_488_858_380_731_577_417_728_000_i128);
        let c33 = - BigInt::from(441_206_965_512_914_835_246_100_i128);
        assert_eq_str!(c33, "-441206965512914835246100");
        let c40 = BigInt::from(1_284_733_132_841_424_456_253_440_i128);
        let c41 = BigInt::from(128_541_798_906_828_816_384_000_i128);
        let c42 = BigInt::from(383_083_609_779_811_215_375_i128);
        let pol = 
          TermBuilder::new().coef(1).xpow(6).build()
        + TermBuilder::new().coef(1).ypow(6).build()
        + TermBuilder::new().coef(-1).xpow(5).ypow(5).build()
        + TermBuilder::new().coef(3_720).xpow(5).ypow(4).build()
        + TermBuilder::new().coef(3_720).ypow(5).xpow(4).build()
        + TermBuilder::new().coef(-4_550_940).xpow(5).ypow(3).build()
        + TermBuilder::new().coef(-4_550_940).ypow(5).xpow(3).build()
        + TermBuilder::new().coef(2_028_551_200).xpow(5).ypow(2).build()
        + TermBuilder::new().coef(2_028_551_200).ypow(5).xpow(2).build()
        + TermBuilder::new().coef(-246_683_410_950).xpow(5).ypow(1).build()
        + TermBuilder::new().coef(-246_683_410_950).xpow(1).ypow(5).build()
        + TermBuilder::new().coef(1_963_211_489_280).xpow(5).build()
        + TermBuilder::new().coef(1_963_211_489_280).ypow(5).build()
        + TermBuilder::new().coef(1_665_999_364_600).xpow(4).ypow(4).build()
        + TermBuilder::new().coef(107_878_928_185_336_800).xpow(4).ypow(3).build()
        + TermBuilder::new().coef(107_878_928_185_336_800).xpow(3).ypow(4).build()
        + TermBuilder::new().coef(&c42).xpow(4).ypow(2).build()
        + TermBuilder::new().coef(&c42).xpow(2).ypow(4).build()
        + TermBuilder::new().coef(&c41).xpow(4).ypow(1).build()
        + TermBuilder::new().coef(&c41).xpow(1).ypow(4).build()
        + TermBuilder::new().coef(&c40).xpow(4).build()
        + TermBuilder::new().coef(&c40).ypow(4).build()
        + TermBuilder::new().coef(&c33).xpow(3).ypow(3).build()
        + TermBuilder::new().coef(&c32).xpow(3).ypow(2).build()
        + TermBuilder::new().coef(&c32).ypow(3).xpow(2).build()
        + TermBuilder::new().coef(&c31).xpow(3).ypow(1).build()
        + TermBuilder::new().coef(&c31).xpow(1).ypow(3).build()
        + TermBuilder::new().coef(&c30).xpow(3).build()
        + TermBuilder::new().coef(&c30).ypow(3).build()
        + TermBuilder::new().coef(&c22).xpow(2).ypow(2).build()
        + TermBuilder::new().coef(&c21).xpow(2).ypow(1).build()
        + TermBuilder::new().coef(&c21).ypow(2).xpow(1).build()
        + TermBuilder::new().coef(&c20).xpow(2).build()
        + TermBuilder::new().coef(&c20).ypow(2).build()
        + TermBuilder::new().coef(&c11).xpow(1).ypow(1).build()
        + TermBuilder::new().coef(&c10).xpow(1).build()
        + TermBuilder::new().coef(&c10).ypow(1).build()
        + TermBuilder::new().coef(&c00).build();
        return pol;
    } else {
        panic!();
    }
}

/// calculate modular polynomial
pub fn subscripted_variable_modular_polynomial(p: u64) -> Polynomial {
    assert!(p >= 2);
    assert!(primes::is_prime(p));
    let mut pol = term_builder::TermBuilder::new()
                .coef(&BigInt::from(1))
                .xpow(&BigInt::from(p+1))
                .build().to_pol();
    pol += term_builder::TermBuilder::new()
                .coef(&BigInt::from(1))
                .ypow(&BigInt::from(p+1))
                .build().to_pol();
    pol += term_builder::TermBuilder::new()
            .coef(&BigInt::from(-1))
            .xpow(&BigInt::from(p))
            .ypow(&BigInt::from(p))
            .build().to_pol();
    pol += term_builder::TermBuilder::new()
            .coef(&BigInt::from(-1))
            .xpow(&BigInt::from(1))
            .ypow(&BigInt::from(1))
            .build().to_pol();

    for i in num_iter::range(0, p) {
        pol += term_builder::TermBuilder::new()
            .coef(&BigInt::from(1))
            .xpow(&BigInt::from(i))
            .ypow(&BigInt::from(i))
            .variable_ij(i, i)
            .build().to_pol();
    }
    for i in num_iter::range(0, p+1) {
        for j in num_iter::range(i+1, p+1) {
            pol += term_builder::TermBuilder::new()
                .coef(&BigInt::from(1))
                .xpow(&BigInt::from(i))
                .ypow(&BigInt::from(j))
                .variable_ij(i, j)
                .build().to_pol();
            pol += term_builder::TermBuilder::new()
                .coef(&BigInt::from(1))
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
    let min = BigInt::from(- pp.power(2) - pp + BigInt::from(1));
    let mut v: Vec<polynomial::Polynomial> = Vec::new();
    for i in num_iter::range(min, BigInt::from(1)) {
        v.push(pol_q.to_q_power_coef(&i));
    }
    v
}

pub fn subscripted_variable_modular_polynomial_q(p: u64) -> polynomial::Polynomial {
    let pol = subscripted_variable_modular_polynomial(p);
    // j-invariant q order needs p^2 + p
    let p = BigInt::from(p);
    let j_order = p.power(2) + p.clone();
    let j = j_invariant::j_invariant(&j_order);
    let j_q_power = j.to_q_power(&p.clone());
    let pol1 = pol.eval_x_polynomial(&j);
    let pol2 = pol1.eval_y_polynomial(&j_q_power);
    let pol_q = pol2.omit_high_order_q(&BigInt::from(0)); 
    pol_q
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
    assert_eq_str!(pol, "x^3 - x^2 y^2 + x^2 y c_1_2 + x^2 c_0_2 + x y^2 c_1_2 + x y c_1_1 - x y + x c_0_1 + y^3 + y^2 c_0_2 + y c_0_1 + c_0_0");

    let pol_q = subscripted_variable_modular_polynomial_q(p);
    assert_eq_str!(pol_q, "126112980648 c_1_2 + 22047296 c_1_1 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 - 941844902299520 + 1495268650 q^-1 c_1_2 + 197628 q^-1 c_1_1 + 1488 q^-1 c_0_2 + q^-1 c_0_1 - 10291427447328 q^-1 + 23548880 q^-2 c_1_2 + 744 q^-2 c_1_1 + 1489 q^-2 c_0_2 + q^-2 c_0_1 - 73882907184 q^-2 + 199860 q^-3 c_1_2 + q^-3 c_1_1 - 338165056 q^-3 + 745 q^-4 c_1_2 + q^-4 c_0_2 - 946560 q^-4 + q^-5 c_1_2 - 1488 q^-5");

    let list = subscripted_variable_modular_polynomial_list(p);
    assert_eq!(list.len(), 6);
    assert_eq_str!(list[0], "c_1_2 - 1488");
    assert_eq_str!(list[1], "745 c_1_2 + c_0_2 - 946560");
    assert_eq_str!(list[2], "199860 c_1_2 + c_1_1 - 338165056");
    assert_eq_str!(list[3], "23548880 c_1_2 + 744 c_1_1 + 1489 c_0_2 + c_0_1 - 73882907184");
    assert_eq_str!(list[4], "1495268650 c_1_2 + 197628 c_1_1 + 1488 c_0_2 + c_0_1 - 10291427447328");
    assert_eq_str!(list[5], "126112980648 c_1_2 + 22047296 c_1_1 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 - 941844902299520");
}

#[test]
fn subscripted_variable_modular_polynomial_p3_test() {
    let p = 3;
    let pol = subscripted_variable_modular_polynomial(p);
    assert_eq_str!(pol, "x^4 - x^3 y^3 + x^3 y^2 c_2_3 + x^3 y c_1_3 + x^3 c_0_3 + x^2 y^3 c_2_3 + x^2 y^2 c_2_2 + x^2 y c_1_2 + x^2 c_0_2 + x y^3 c_1_3 + x y^2 c_1_2 + x y c_1_1 - x y + x c_0_1 + y^4 + y^3 c_0_3 + y^2 c_0_2 + y c_0_1 + c_0_0");

    let pol_q = subscripted_variable_modular_polynomial_q(p);
    assert_eq_str!(pol_q, "65201219062152355627872 c_2_3 + 1769523054225123330 c_2_2 + 45070853378831784 c_1_3 + 15329636199360 c_1_2 + 864853506 c_1_1 + 2710404480 c_0_3 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 - 4197391688509447725659608686 + 1768996792073642463801 q^-1 c_2_3 + 45738435798911040 q^-1 c_2_2 + 1345112306562240 q^-1 c_1_3 + 437662034132 q^-1 c_1_2 + 21494504 q^-1 c_1_1 + 2251260 q^-1 c_0_3 + 1488 q^-1 c_0_2 + q^-1 c_0_1 - 109642480943384827151178336 q^-1 + 57383938307316236432 q^-2 c_2_3 + 833607524819048 q^-2 c_2_2 + 90830762088165 q^-2 c_1_3 + 20874771304 q^-2 c_1_2 + 196884 q^-2 c_1_1 + 2232 q^-2 c_0_3 + q^-2 c_0_2 - 2319314637883533600896496 q^-2 + 1807266275628309954 q^-3 c_2_3 + 9791534543904 q^-3 c_2_2 + 6184170973560 q^-3 c_1_3 + 866354346 q^-3 c_1_2 + 744 q^-3 c_1_1 + 2251261 q^-3 c_0_3 + 1488 q^-3 c_0_2 + q^-3 c_0_1 - 38360228325192295302456 q^-3 + 45943884263343552 q^-4 c_2_3 + 72476838420 q^-4 c_2_2 + 381181215440 q^-4 c_1_3 + 21496736 q^-4 c_1_2 + q^-4 c_1_1 - 480257921200323487545 q^-4 + 834399255041138 q^-5 c_2_3 + 335952400 q^-5 c_2_2 + 20685303576 q^-5 c_1_3 + 196885 q^-5 c_1_2 - 4461984465203723760 q^-5 + 9793594541808 q^-6 c_2_3 + 947304 q^-6 c_2_2 + 865960579 q^-6 c_1_3 + 744 q^-6 c_1_2 + 2232 q^-6 c_0_3 + q^-6 c_0_2 - 30459141461183076 q^-6 + 72480196752 q^-7 c_2_3 + 1488 q^-7 c_2_2 + 21495992 q^-7 c_1_3 + q^-7 c_1_2 - 151527078622080 q^-7 + 335955376 q^-8 c_2_3 + q^-8 c_2_2 + 196884 q^-8 c_1_3 - 541783100214 q^-8 + 947305 q^-9 c_2_3 + 744 q^-9 c_1_3 + q^-9 c_0_3 - 1355201496 q^-9 + 1488 q^-10 c_2_3 + q^-10 c_1_3 - 2251260 q^-10 + q^-11 c_2_3 - 2232 q^-11");

    let list = subscripted_variable_modular_polynomial_list(p);
    assert_eq!(list.len(), 12);
    assert_eq_str!(list[0], "c_2_3 - 2232");
    assert_eq_str!(list[1], "1488 c_2_3 + c_1_3 - 2251260");
    assert_eq_str!(list[2], "947305 c_2_3 + 744 c_1_3 + c_0_3 - 1355201496");
    assert_eq_str!(list[3], "335955376 c_2_3 + c_2_2 + 196884 c_1_3 - 541783100214");
    assert_eq_str!(list[4], "72480196752 c_2_3 + 1488 c_2_2 + 21495992 c_1_3 + c_1_2 - 151527078622080");
    assert_eq_str!(list[5], "9793594541808 c_2_3 + 947304 c_2_2 + 865960579 c_1_3 + 744 c_1_2 + 2232 c_0_3 + c_0_2 - 30459141461183076");
    assert_eq_str!(list[6], "834399255041138 c_2_3 + 335952400 c_2_2 + 20685303576 c_1_3 + 196885 c_1_2 - 4461984465203723760");
    assert_eq_str!(list[7], "45943884263343552 c_2_3 + 72476838420 c_2_2 + 381181215440 c_1_3 + 21496736 c_1_2 + c_1_1 - 480257921200323487545");
    assert_eq_str!(list[8], "1807266275628309954 c_2_3 + 9791534543904 c_2_2 + 6184170973560 c_1_3 + 866354346 c_1_2 + 744 c_1_1 + 2251261 c_0_3 + 1488 c_0_2 + c_0_1 - 38360228325192295302456");
    assert_eq_str!(list[9], "57383938307316236432 c_2_3 + 833607524819048 c_2_2 + 90830762088165 c_1_3 + 20874771304 c_1_2 + 196884 c_1_1 + 2232 c_0_3 + c_0_2 - 2319314637883533600896496");
    assert_eq_str!(list[10], "1768996792073642463801 c_2_3 + 45738435798911040 c_2_2 + 1345112306562240 c_1_3 + 437662034132 c_1_2 + 21494504 c_1_1 + 2251260 c_0_3 + 1488 c_0_2 + c_0_1 - 109642480943384827151178336");
    assert_eq_str!(list[11], "65201219062152355627872 c_2_3 + 1769523054225123330 c_2_2 + 45070853378831784 c_1_3 + 15329636199360 c_1_2 + 864853506 c_1_1 + 2710404480 c_0_3 + 1894608 c_0_2 + 1488 c_0_1 + c_0_0 - 4197391688509447725659608686");
}

