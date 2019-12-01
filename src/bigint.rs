use num_integer::Integer;
use num_bigint::BigInt;
use num_traits::Zero;
use num_traits::One;

/// T^n
/// NOTE: BigInt::Pow is not enough functionality, so implement by myself.
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

impl Power<BigInt> for BigInt { 
    fn power(&self, n: BigInt) -> Self {
        self.power(&n)
    }
}

impl<'a> Power<&'a BigInt> for BigInt { 
    fn power(&self, n: &BigInt) -> Self {
        if n == &Zero::zero() {
            return One::one();
        }
        let mut e = n.clone();
        let mut b = self.clone();
        let mut r = BigInt::from(1);
        while e > One::one() {
            if e.is_odd() {
                r *= b.clone();
            }
            let f = b.clone() * b.clone();
            b = f;
            e /= 2;
        } 
        r * b
    }
}

impl Power<i64> for BigInt {
    fn power(&self, n: i64) -> Self {
        if n == 0 {
            return One::one();
        }
        let mut e = n;
        let mut b = self.clone();
        let mut r = BigInt::from(1);
        while e > One::one() {
            if e.is_odd() {
                r *= b.clone();
            }
            let f = b.clone() * b.clone();
            b = f;
            e /= 2;
        }
        r * b
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

/// x (mod l) = r
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ModResult {
    pub l: BigInt,
    pub r: BigInt,
}

/// Chinese remainder thereom
/// if x mod l1 = r1, x mod l2 = r2, ... x mod ln = rn
/// then x mod (l1 * l2 * ... * ln) = R
///
/// return: l1 * l2 * ... * ln, R
pub fn chinese_remainder(mod_result: &Vec<ModResult>) -> ModResult {
    let mut l: BigInt = One::one();
    let mut r: BigInt = One::one();
    for res in mod_result.iter() {
        let (_, p, _) = extended_gcd(l.clone(), res.l.clone());
        r += (res.r.clone() - r.clone()) * l.clone() * p;
        l *= res.l.clone();
        println!("r:{} l:{}", r, l);
    }
    ModResult { l, r }
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
    assert_eq!(Integer::lcm(&BigInt::from(4), &BigInt::from(6)), BigInt::from(12));
    assert_eq!(Integer::lcm(&BigInt::from(-4), &BigInt::from(6)), BigInt::from(12));
    assert_eq!(Integer::lcm(&BigInt::from(-4), &BigInt::from(-6)), BigInt::from(12));
}

#[test]
fn chinese_remainder_test() {
    use self::chinese_remainder;
    let mod_result = vec![
        ModResult {
            l: BigInt::from(3),
            r: BigInt::from(2),
        },
        ModResult {
            l: BigInt::from(5),
            r: BigInt::from(3),
        },
        ModResult {
            l: BigInt::from(7),
            r: BigInt::from(2),
        }];
    assert_eq!(chinese_remainder(&mod_result), 
               ModResult { 
                   l: BigInt::from(3*5*7), 
                   r: BigInt::from(23 - 3*5*7) 
               });
}

