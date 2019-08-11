use num_integer::Integer;
use num_bigint::BigInt;


pub trait Power<T> {
    fn power(&self, n: T) -> Self;
}

pub trait PowerModular {
    fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self;
}

pub trait Inverse {
    fn inverse(&self, p: &BigInt) -> Self;
}

impl<'a> Power<&'a BigInt> for BigInt { 
    fn power(&self, n: &BigInt) -> Self {
        let mut t = BigInt::from(1);
        for _i in num_iter::range(BigInt::from(0), n.clone()) {
            t *= self;
        }
        t
    }
}

impl Power<i64> for BigInt {
    fn power(&self, n: i64) -> Self {
        let n = BigInt::from(n);
        self.power(&n)
    }
}

impl PowerModular for BigInt {
    fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self {
        let mut t = BigInt::from(1);
        for _i in num_iter::range(BigInt::from(0), n.clone()) {
            t *= self;
        }
        t = t.mod_floor(p);
        t
    }
}

impl Inverse for BigInt {
    fn inverse(&self, p: &BigInt) -> Self {
        assert!(p >= &BigInt::from(2));
        self.power_modular(&(p.clone()-&BigInt::from(2)), &p)
    }
}

#[test]
fn bigint_power_test() {
    let q = BigInt::from(2);
    let q2 = q.power(&BigInt::from(3));
    assert_eq_str!(q2, "8");
}

#[test]
fn bigint_power_modular_test() {
    let p = BigInt::from(5);
    assert_eq_str!(BigInt::from(3).power_modular(&BigInt::from(2), &p), "4");
    assert_eq_str!(BigInt::from(3).power_modular(&BigInt::from(3), &p), "2");
    assert_eq_str!(BigInt::from(3).power_modular(&BigInt::from(4), &p), "1");
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

