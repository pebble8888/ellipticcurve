extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;
extern crate divrem;

use num_bigint::BigInt;
use num_traits::Zero;

pub trait Power {
    fn power(&self, n: &BigInt) -> Self; 
    fn power_i(&self, n: i64) -> Self;
}

pub trait PowerModular {
    fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self;
}

pub trait Inverse {
    fn inverse(&self, p: &BigInt) -> Self;
}

pub trait DivFloor<Rhs=Self> {
    type Output;
    fn div_floor(&self, rhs: &Rhs) -> Self::Output;
}

pub trait RemFloor<Rhs=Self> {
    type Output;
    fn rem_floor(&self, rhs: &Rhs) -> Self::Output;
}

impl Power for BigInt { 
    fn power(&self, n: &BigInt) -> Self {
        let mut t = BigInt::from(1);
        for _i in num_iter::range(BigInt::from(0), n.clone()) {
            t = &t * self;
        }
        t
    }

    fn power_i(&self, n: i64) -> Self {
        let n = BigInt::from(n);
        self.power(&n)
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

impl DivFloor for BigInt {
    type Output = Self;
    fn div_floor(&self, other: &Self) -> Self {
        if self > &Zero::zero() && other > &Zero::zero() {
            self / other
        } else if self > &Zero::zero() && other < &Zero::zero() {
            ((self - BigInt::from(1)) / other) - BigInt::from(1)
        } else if self < &Zero::zero() && other > &Zero::zero() {
            ((self + BigInt::from(1)) / other) - BigInt::from(1)
        } else {
            self / other
        }
    }
}

impl RemFloor for BigInt {
    type Output = Self;
    fn rem_floor(&self, other: &Self) -> Self {
        if self > &Zero::zero() && other > &Zero::zero() {
            self % other
        } else if self > &Zero::zero() && other < &Zero::zero() {
            ((self - BigInt::from(1)) % other) + other + BigInt::from(1)
        } else if self < &Zero::zero() && other > &Zero::zero() {
            ((self + BigInt::from(1)) % other) + other - BigInt::from(1)
        } else {
            self % other
        }
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
    let p = BigInt::from(19);
    assert_eq!(BigInt::from(1).inverse(&p).to_string(), "1");
    assert_eq!(BigInt::from(2).inverse(&p).to_string(), "10");
    assert_eq!(BigInt::from(3).inverse(&p).to_string(), "13");
    assert_eq!(BigInt::from(4).inverse(&p).to_string(), "5");
}

//use divrem::DivFloor;
//use divrem::RemFloor;

#[test]
fn bigint_divide_test() {
    /*
    assert_eq!((BigInt::from(7) / BigInt::from(2)).to_string(), "3"); 
    assert_eq!((BigInt::from(7) / BigInt::from(-2)).to_string(), "-3"); 
    assert_eq!((BigInt::from(-7) / BigInt::from(2)).to_string(), "-3"); 
    assert_eq!((BigInt::from(-7) / BigInt::from(-2)).to_string(), "3"); 

    assert_eq!((BigInt::from(7) % BigInt::from(3)).to_string(), "1"); 
    assert_eq!((BigInt::from(7) % BigInt::from(-3)).to_string(), "1"); 
    assert_eq!((BigInt::from(-7) % BigInt::from(3)).to_string(), "-1"); 
    assert_eq!((BigInt::from(-7) % BigInt::from(-3)).to_string(), "-1"); 

    assert_eq!(7.div_floor(2), 3);
    assert_eq!(7.div_floor(-2), -4);
    assert_eq!((-7).div_floor(2), -4);
    assert_eq!((-7).div_floor(-2), 3);

    assert_eq!(7.rem_floor(3), 1);
    assert_eq!(7.rem_floor(-3), -2);
    assert_eq!((-7).rem_floor(3), 2);
    assert_eq!((-7).rem_floor(-3), -1);
    */

    assert_eq!((BigInt::from(7).div_floor(&BigInt::from(2))).to_string(), "3"); 
    assert_eq!((BigInt::from(7).div_floor(&BigInt::from(-2))).to_string(), "-4"); 
    assert_eq!((BigInt::from(-7).div_floor(&BigInt::from(2))).to_string(), "-4"); 
    assert_eq!((BigInt::from(-7).div_floor(&BigInt::from(-2))).to_string(), "3"); 

    assert_eq!((BigInt::from(7).rem_floor(&BigInt::from(3))).to_string(), "1"); 
    assert_eq!((BigInt::from(7).rem_floor(&BigInt::from(-3))).to_string(), "-2"); 
    assert_eq!((BigInt::from(-7).rem_floor(&BigInt::from(3))).to_string(), "2"); 
    assert_eq!((BigInt::from(-7).rem_floor(&BigInt::from(-3))).to_string(), "-1"); 
}

