use num_integer::Integer;
use num_bigint::BigInt;
use num_traits::Zero;
use num_traits::One;

/// T^n
pub trait Power<T> {
    fn power(&self, n: T) -> Self;
}

/// T^n (mod p)
pub trait PowerModulo {
    fn power_modulo(&self, n: &BigInt, p: &BigInt) -> Self;
}

/// 1/T (mod p)
pub trait Inverse {
    fn inverse(&self, p: &BigInt) -> Self;
}

impl<'a> Power<&'a BigInt> for BigInt { 
    fn power(&self, n: &BigInt) -> Self {
        let mut t: BigInt = One::one();
        for _i in num_iter::range(Zero::zero(), n.clone()) {
            t *= self;
        }
        t
    }
}

impl Power<i64> for BigInt {
    fn power(&self, n: i64) -> Self {
        let mut t: BigInt = One::one();
        for _i in num_iter::range(Zero::zero(), n.clone()) {
            t *= self;
        }
        t
    }
}

impl PowerModulo for BigInt {
    fn power_modulo(&self, n: &BigInt, p: &BigInt) -> Self {
        let mut t = BigInt::from(1);
        for _i in num_iter::range(Zero::zero(), n.clone()) {
            t *= self;
        }
        t = t.mod_floor(p);
        t
    }
}

impl Inverse for BigInt {
    fn inverse(&self, p: &BigInt) -> Self {
        assert!(p >= &BigInt::from(2));
        self.power_modulo(&(p.clone()-&BigInt::from(2)), &p)
    }
}

/// extended euclid algorithm
/// 
/// a x + b y = gcd(a, b)
///
/// params: (a, b) 
/// return: (gcd(a, b), x, y)
///
pub fn extended_gcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if a.clone() == Zero::zero() {
        return (b, Zero::zero(), One::one());
    }
    let r = b.clone() % a.clone();
    let q = b.clone() / a.clone();
    let (g, x, y) = extended_gcd(r, a);
    (g, y - q * x.clone(), x.clone())
}



#[test]
fn bigint_power_test() {
    let q = BigInt::from(2);
    let q2 = q.power(&BigInt::from(3));
    assert_eq_str!(q2, "8");
}

#[test]
fn bigint_power_modulo_test() {
    let p = BigInt::from(5);
    assert_eq_str!(BigInt::from(3).power_modulo(&BigInt::from(2), &p), "4");
    assert_eq_str!(BigInt::from(3).power_modulo(&BigInt::from(3), &p), "2");
    assert_eq_str!(BigInt::from(3).power_modulo(&BigInt::from(4), &p), "1");
}

#[test]
fn bigint_inverse_test() {
    let p = BigInt::from(19);
    assert_eq_str!(BigInt::from(1).inverse(&p), "1");
    assert_eq_str!(BigInt::from(2).inverse(&p), "10");
    assert_eq_str!(BigInt::from(3).inverse(&p), "13");
    assert_eq_str!(BigInt::from(4).inverse(&p), "5");
}

#[test]
fn bigint_divide_test() {
    assert_eq!((BigInt::from(7) / BigInt::from(2)).to_string(), "3"); 
    assert_eq!((BigInt::from(7) / BigInt::from(-2)).to_string(), "-3"); 
    assert_eq!((BigInt::from(-7) / BigInt::from(2)).to_string(), "-3"); 
    assert_eq!((BigInt::from(-7) / BigInt::from(-2)).to_string(), "3"); 

    assert_eq!((BigInt::from(7) % BigInt::from(3)).to_string(), "1"); 
    assert_eq!((BigInt::from(7) % BigInt::from(-3)).to_string(), "1"); 
    assert_eq!((BigInt::from(-7) % BigInt::from(3)).to_string(), "-1"); 
    assert_eq!((BigInt::from(-7) % BigInt::from(-3)).to_string(), "-1"); 

    assert_eq!(7.div_floor(&2), 3);
    assert_eq!(7.div_floor(&-2), -4);
    assert_eq!((-7).div_floor(&2), -4);
    assert_eq!((-7).div_floor(&-2), 3);

    assert_eq!(7.mod_floor(&3), 1);
    assert_eq!(7.mod_floor(&-3), -2);
    assert_eq!((-7).mod_floor(&3), 2);
    assert_eq!((-7).mod_floor(&-3), -1);

    assert_eq_str!((BigInt::from(7).div_floor(&BigInt::from(2))), "3"); 
    assert_eq_str!((BigInt::from(7).div_floor(&BigInt::from(-2))), "-4"); 
    assert_eq_str!((BigInt::from(-7).div_floor(&BigInt::from(2))), "-4"); 
    assert_eq_str!((BigInt::from(-7).div_floor(&BigInt::from(-2))), "3"); 

    assert_eq_str!((BigInt::from(7).mod_floor(&BigInt::from(3))), "1"); 
    assert_eq_str!((BigInt::from(7).mod_floor(&BigInt::from(-3))), "-2"); 
    assert_eq_str!((BigInt::from(-7).mod_floor(&BigInt::from(3))), "2"); 
    assert_eq_str!((BigInt::from(-7).mod_floor(&BigInt::from(-3))), "-1"); 
}

#[test]
fn extended_euclid_test1() {
    let (g, x, y) = extended_gcd(BigInt::from(8), BigInt::from(11)); 
    assert_eq!(g, BigInt::from(1));
    assert_eq!(x, BigInt::from(-4));
    assert_eq!(y, BigInt::from(3));
}

#[test]
fn extended_euclid_test2() {
    let (g, x, y) = extended_gcd(BigInt::from(11), BigInt::from(8)); 
    assert_eq!(g, BigInt::from(1));
    assert_eq!(x, BigInt::from(3));
    assert_eq!(y, BigInt::from(-4));
}

#[test]
fn lcm_test() {
    let a = BigInt::from(4);
    assert_eq!(Integer::lcm(&BigInt::from(4), &BigInt::from(6)), BigInt::from(12));
    assert_eq!(Integer::lcm(&BigInt::from(-4), &BigInt::from(6)), BigInt::from(12));
    assert_eq!(Integer::lcm(&BigInt::from(-4), &BigInt::from(-6)), BigInt::from(12));
}

