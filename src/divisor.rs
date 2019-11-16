use num_bigint::BigInt;
use num_traits::Zero;
use crate::bigint::Power;

/// The sum of positive divisors function 
pub fn sigma_divisor(n: &BigInt, power: &BigInt) -> BigInt {
    if *power < Zero::zero() {
        return Zero::zero();
    }
    if *n < Zero::zero() {
        return Zero::zero();
    }
    let mut sum = BigInt::from(0);
    let n_plus_1 = n.clone() + BigInt::from(1);
    for i in num_iter::range(BigInt::from(1), n_plus_1) {
        if n % &i == Zero::zero() {
            sum += i.power(power);
        }
    }
    sum
}

#[test]
fn sigma_divisor_test1() {
    assert_eq_str!(sigma_divisor(&BigInt::from(1), &BigInt::from(0)), "1");
    assert_eq_str!(sigma_divisor(&BigInt::from(2), &BigInt::from(0)), "2");
    assert_eq_str!(sigma_divisor(&BigInt::from(3), &BigInt::from(0)), "2");
    assert_eq_str!(sigma_divisor(&BigInt::from(4), &BigInt::from(0)), "3");
    assert_eq_str!(sigma_divisor(&BigInt::from(5), &BigInt::from(0)), "2");
    assert_eq_str!(sigma_divisor(&BigInt::from(6), &BigInt::from(0)), "4");
    assert_eq_str!(sigma_divisor(&BigInt::from(7), &BigInt::from(0)), "2");
    assert_eq_str!(sigma_divisor(&BigInt::from(8), &BigInt::from(0)), "4");
    assert_eq_str!(sigma_divisor(&BigInt::from(9), &BigInt::from(0)), "3");
    assert_eq_str!(sigma_divisor(&BigInt::from(10), &BigInt::from(0)), "4");
    assert_eq_str!(sigma_divisor(&BigInt::from(11), &BigInt::from(0)), "2");
    assert_eq_str!(sigma_divisor(&BigInt::from(12), &BigInt::from(0)), "6");
    assert_eq_str!(sigma_divisor(&BigInt::from(13), &BigInt::from(0)), "2");

    assert_eq_str!(sigma_divisor(&BigInt::from(1), &BigInt::from(1)), "1");
    assert_eq_str!(sigma_divisor(&BigInt::from(2), &BigInt::from(1)), "3");
    assert_eq_str!(sigma_divisor(&BigInt::from(3), &BigInt::from(1)), "4");
    assert_eq_str!(sigma_divisor(&BigInt::from(4), &BigInt::from(1)), "7");
    assert_eq_str!(sigma_divisor(&BigInt::from(5), &BigInt::from(1)), "6");
    assert_eq_str!(sigma_divisor(&BigInt::from(6), &BigInt::from(1)), "12");
    assert_eq_str!(sigma_divisor(&BigInt::from(7), &BigInt::from(1)), "8");
    assert_eq_str!(sigma_divisor(&BigInt::from(8), &BigInt::from(1)), "15");
    assert_eq_str!(sigma_divisor(&BigInt::from(9), &BigInt::from(1)), "13");
    assert_eq_str!(sigma_divisor(&BigInt::from(10), &BigInt::from(1)), "18");
    assert_eq_str!(sigma_divisor(&BigInt::from(11), &BigInt::from(1)), "12");
    assert_eq_str!(sigma_divisor(&BigInt::from(12), &BigInt::from(1)), "28");
    assert_eq_str!(sigma_divisor(&BigInt::from(13), &BigInt::from(1)), "14");

    assert_eq_str!(sigma_divisor(&BigInt::from(1), &BigInt::from(3)), "1");
    assert_eq_str!(sigma_divisor(&BigInt::from(2), &BigInt::from(3)), "9");
    assert_eq_str!(sigma_divisor(&BigInt::from(3), &BigInt::from(3)), "28");
    assert_eq_str!(sigma_divisor(&BigInt::from(4), &BigInt::from(3)), "73");
    assert_eq_str!(sigma_divisor(&BigInt::from(5), &BigInt::from(3)), "126");
    assert_eq_str!(sigma_divisor(&BigInt::from(6), &BigInt::from(3)), "252");
    assert_eq_str!(sigma_divisor(&BigInt::from(7), &BigInt::from(3)), "344");
    assert_eq_str!(sigma_divisor(&BigInt::from(8), &BigInt::from(3)), "585");
    assert_eq_str!(sigma_divisor(&BigInt::from(9), &BigInt::from(3)), "757");
    assert_eq_str!(sigma_divisor(&BigInt::from(10), &BigInt::from(3)), "1134");
    assert_eq_str!(sigma_divisor(&BigInt::from(11), &BigInt::from(3)), "1332");
    assert_eq_str!(sigma_divisor(&BigInt::from(12), &BigInt::from(3)), "2044");
    assert_eq_str!(sigma_divisor(&BigInt::from(13), &BigInt::from(3)), "2198");
}

