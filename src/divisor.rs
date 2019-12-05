use num_bigint::BigInt;
use num_traits::{Zero};
use crate::bigint::Power;

/// The sum of positive divisors function 
/// sigma_{power}(n)
pub fn sigma_divisor(n: i32, power: i32) -> BigInt {
    if power < Zero::zero() {
        return Zero::zero();
    }
    if n < Zero::zero() {
        return Zero::zero();
    }
    let mut sum: BigInt = Zero::zero();
    for i in num_iter::range(1, n + 1) {
        if n % i == Zero::zero() {
            sum += BigInt::from(i).power(power);
        }
    }
    sum
}

#[test]
fn sigma_divisor_test1() {
    assert_eq_str!(sigma_divisor(1, 0), "1");
    assert_eq_str!(sigma_divisor(2, 0), "2");
    assert_eq_str!(sigma_divisor(3, 0), "2");
    assert_eq_str!(sigma_divisor(4, 0), "3");
    assert_eq_str!(sigma_divisor(5, 0), "2");
    assert_eq_str!(sigma_divisor(6, 0), "4");
    assert_eq_str!(sigma_divisor(7, 0), "2");
    assert_eq_str!(sigma_divisor(8, 0), "4");
    assert_eq_str!(sigma_divisor(9, 0), "3");
    assert_eq_str!(sigma_divisor(10, 0), "4");
    assert_eq_str!(sigma_divisor(11, 0), "2");
    assert_eq_str!(sigma_divisor(12, 0), "6");
    assert_eq_str!(sigma_divisor(13, 0), "2");

    assert_eq_str!(sigma_divisor(1, 1), "1");
    assert_eq_str!(sigma_divisor(2, 1), "3");
    assert_eq_str!(sigma_divisor(3, 1), "4");
    assert_eq_str!(sigma_divisor(4, 1), "7");
    assert_eq_str!(sigma_divisor(5, 1), "6");
    assert_eq_str!(sigma_divisor(6, 1), "12");
    assert_eq_str!(sigma_divisor(7, 1), "8");
    assert_eq_str!(sigma_divisor(8, 1), "15");
    assert_eq_str!(sigma_divisor(9, 1), "13");
    assert_eq_str!(sigma_divisor(10, 1), "18");
    assert_eq_str!(sigma_divisor(11, 1), "12");
    assert_eq_str!(sigma_divisor(12, 1), "28");
    assert_eq_str!(sigma_divisor(13, 1), "14");

    assert_eq_str!(sigma_divisor(1, 3), "1");
    assert_eq_str!(sigma_divisor(2, 3), "9");
    assert_eq_str!(sigma_divisor(3, 3), "28");
    assert_eq_str!(sigma_divisor(4, 3), "73");
    assert_eq_str!(sigma_divisor(5, 3), "126");
    assert_eq_str!(sigma_divisor(6, 3), "252");
    assert_eq_str!(sigma_divisor(7, 3), "344");
    assert_eq_str!(sigma_divisor(8, 3), "585");
    assert_eq_str!(sigma_divisor(9, 3), "757");
    assert_eq_str!(sigma_divisor(10, 3), "1134");
    assert_eq_str!(sigma_divisor(11, 3), "1332");
    assert_eq_str!(sigma_divisor(12, 3), "2044");
    assert_eq_str!(sigma_divisor(13, 3), "2198");
}

