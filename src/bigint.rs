extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;

pub trait Power {
    fn power(&self, n: &BigInt) -> Self; 
}

pub trait PowerModular {
    fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self;
}

pub trait Inverse {
    fn inverse(&self, p: &BigInt) -> Self;
}

impl Power for BigInt { 
    fn power(&self, n: &BigInt) -> Self {
        let mut t = BigInt::from(1);
        for _i in num_iter::range(BigInt::from(0), n.clone()) {
            t = &t * self;
        }
        t
    }
}

impl PowerModular for BigInt {
    fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self {
        let mut t = BigInt::from(1);
        for _i in num_iter::range(BigInt::from(0), n.clone()) {
            t = &t * self;
            t = &t % p;
        }
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
    assert_eq!(q2.to_string(), "8");
}

#[test]
fn bigint_power_modular_test() {
    let p = BigInt::from(5);
    assert_eq!(BigInt::from(3).power_modular(&BigInt::from(2), &p).to_string(), "4");
    assert_eq!(BigInt::from(3).power_modular(&BigInt::from(3), &p).to_string(), "2");
    assert_eq!(BigInt::from(3).power_modular(&BigInt::from(4), &p).to_string(), "1");
}

#[test]
fn bigint_inverse_test() {
    let p = BigInt::from(5);
    assert_eq!(BigInt::from(1).inverse(&p).to_string(), "1");
    assert_eq!(BigInt::from(2).inverse(&p).to_string(), "3");
    assert_eq!(BigInt::from(3).inverse(&p).to_string(), "2");
    assert_eq!(BigInt::from(4).inverse(&p).to_string(), "4");
}

