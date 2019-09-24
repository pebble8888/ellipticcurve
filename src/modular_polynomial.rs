use num_bigint::BigInt;
use super::polynomial;
use super::term_builder;
use super::term_builder::TermBuildable;
use super::bigint::Power;

type TermBuilder = term_builder::TermBuilder;
type Polynomial = polynomial::Polynomial;

pub fn modular_polynomial(n: &BigInt) -> Polynomial {
    assert!(*n >= BigInt::from(2));
    if *n == BigInt::from(2) {
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
    } else if *n == BigInt::from(3) {
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
    } else if *n == BigInt::from(5) {
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

#[test]
fn modular_polynomial_test2() {
    let pol = modular_polynomial(&BigInt::from(2));
    assert_eq_str!(pol, "x^3 - x^2 y^2 + 1488 x^2 y - 162000 x^2 + 1488 x y^2 + 40773375 x y + 8748000000 x + y^3 - 162000 y^2 + 8748000000 y - 157464000000000");
}

#[test]
fn modular_polynomial_test3() {
    let pol = modular_polynomial(&BigInt::from(3));
    assert_eq_str!(pol, "x^4 - x^3 y^3 + 2232 x^3 y^2 - 1069956 x^3 y + 36864000 x^3 + 2232 x^2 y^3 + 2587918086 x^2 y^2 + 8900222976000 x^2 y + 452984832000000 x^2 - 1069956 x y^3 + 8900222976000 x y^2 - 770845966336000000 x y + 1855425871872000000000 x + y^4 + 36864000 y^3 + 452984832000000 y^2 + 1855425871872000000000 y");
}

#[test]
fn modular_polynomial_test5() {
    let pol = modular_polynomial(&BigInt::from(5));
    assert_eq_str!(pol, "x^6 - x^5 y^5 + 3720 x^5 y^4 - 4550940 x^5 y^3 + 2028551200 x^5 y^2 - 246683410950 x^5 y + 1963211489280 x^5 + 3720 x^4 y^5 + 1665999364600 x^4 y^4 + 107878928185336800 x^4 y^3 + 383083609779811215375 x^4 y^2 + 128541798906828816384000 x^4 y + 1284733132841424456253440 x^4 - 4550940 x^3 y^5 + 107878928185336800 x^3 y^4 - 441206965512914835246100 x^3 y^3 + 26898488858380731577417728000 x^3 y^2 - 192457934618928299655108231168000 x^3 y + 280244777828439527804321565297868800 x^3 + 2028551200 x^2 y^5 + 383083609779811215375 x^2 y^4 + 26898488858380731577417728000 x^2 y^3 + 5110941777552418083110765199360000 x^2 y^2 + 36554736583949629295706472332656640000 x^2 y + 6692500042627997708487149415015068467200 x^2 - 246683410950 x y^5 + 128541798906828816384000 x y^4 - 192457934618928299655108231168000 x y^3 + 36554736583949629295706472332656640000 x y^2 - 264073457076620596259715790247978782949376 x y + 53274330803424425450420160273356509151232000 x + y^6 + 1963211489280 y^5 + 1284733132841424456253440 y^4 + 280244777828439527804321565297868800 y^3 + 6692500042627997708487149415015068467200 y^2 + 53274330803424425450420160273356509151232000 y + 141359947154721358697753474691071362751004672000");
}
