extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use num_traits::pow;
use std::cmp::Ordering;
use std::ops::Mul;
use std::ops::Div;
use std::ops::Neg;
use std::fmt;
//use std::default::Default;

// coef x^xpow y^ypow
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Unit {
    pub coef: BigInt,
    pub xpow: BigInt,
    pub ypow: BigInt,
}

impl Mul for Unit {
    type Output = Unit;
    fn mul(self, other: Self) -> Self {
        Unit {
            coef: &self.coef * &other.coef,
            xpow: &self.xpow + &other.xpow,
            ypow: &self.ypow + &other.ypow,
        }
    }
}

impl Div for Unit {
    type Output = Unit;
    fn div(self, other: Self) -> Self {
        Unit {
            coef: &self.coef / &other.coef,
            xpow: &self.xpow - &other.xpow,
            ypow: &self.ypow - &other.ypow,
        }
    }
}

impl Neg for Unit {
    type Output = Unit;
    fn neg(self) -> Self {
        Unit {
            coef: -&self.coef,
            xpow: self.xpow,
            ypow: self.ypow,
        }
    }
}

impl Unit {
    pub fn equal_order(&self, other: &Self) -> bool {
        return &self.xpow == &other.xpow && self.ypow == other.ypow 
    }
    pub fn power(&self, val: usize) -> Self {
        let coef = pow(self.coef.clone(), val);
        let xpow = self.xpow.clone() * val;
        let ypow = self.ypow.clone() * val;
        Unit {
            coef: coef,
            xpow: xpow,
            ypow: ypow,
        }
    }
    pub fn to_frob(&self, val: usize) -> Self {
        let coef = self.coef.clone();
        let xpow = self.xpow.clone() * val;
        let ypow = self.ypow.clone() * val;
        Unit {
            coef: coef,
            xpow: xpow,
            ypow: ypow,
        }
    }
    pub fn modular(&self, val: BigInt) -> Self {
        let coef = self.coef.clone() % val;
        Unit {
            coef: coef,
            xpow: self.xpow.clone(),
            ypow: self.ypow.clone(),
        }
    }
    pub fn is_zero(&self) -> bool {
        self.coef == Zero::zero()
    }

    pub fn has_y(&self) -> bool {
        self.ypow != Zero::zero()
    }
}

impl Ord for Unit {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.xpow < other.xpow {
            return Ordering::Less;
        } else if self.xpow > other.xpow {
            return Ordering::Greater;
        }
        if self.ypow < other.ypow {
            return Ordering::Less;
        } else if self.ypow > other.ypow {
            return Ordering::Greater;
        }
        return Ordering::Equal;
    }
}

impl PartialOrd for Unit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Unit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.coef == One::one() {
            if self.xpow == Zero::zero() && self.ypow == Zero::zero() {
                write!(f, "1")
            } else {
                let mut st = String::new();
                if self.xpow != Zero::zero() {
                    if self.xpow == One::one() {
                        st.push_str("x");
                    } else {
                        st.push_str("x^");
                        st.push_str(&self.xpow.to_string());
                    }
                    st.push_str(" ");
                }
                if self.ypow != Zero::zero() {
                    if self.ypow == One::one() {
                        st.push_str("y");
                    } else {
                        st.push_str("y^");
                        st.push_str(&self.ypow.to_string());
                    }
                }
                write!(f, "{}", st.trim_end())
            }
        } else if self.coef == BigInt::from(-1) {
            if self.xpow == Zero::zero() && self.ypow == Zero::zero() {
                write!(f, "- 1")
            } else {
                let mut st = String::new();
                st.push_str("- ");
                if self.xpow != Zero::zero() {
                    if self.xpow == One::one() {
                        st.push_str("x");
                    } else {
                        st.push_str("x^");
                        st.push_str(&self.xpow.to_string());
                    }
                    st.push_str(" ");
                }
                if self.ypow != Zero::zero() {
                    if self.ypow == One::one() {
                        st.push_str("y");
                    } else {
                        st.push_str("y^");
                        st.push_str(&self.ypow.to_string());
                    }
                }
                write!(f, "{}", st.trim_end())
            }

        } else {
            let mut st = String::new();
            if self.coef >= Zero::zero() {
                st.push_str(&self.coef.to_string());
            } else {
                let abs_coef = self.coef.clone() * BigInt::from(-1);
                st.push_str("- ");
                st.push_str(&abs_coef.to_string());
            }
            st.push_str(" ");
            if self.xpow != Zero::zero() {
                if self.xpow == One::one() {
                    st.push_str("x");
                } else {
                    st.push_str("x^");
                    st.push_str(&self.xpow.to_string());
                }
                st.push_str(" ");
            }
            if self.ypow != Zero::zero() {
                if self.ypow == One::one() {
                    st.push_str("y");
                } else {
                    st.push_str("y^");
                    st.push_str(&self.ypow.to_string());
                }
            }
            write!(f, "{}", st.trim_end())
        }
    }
}

#[test]
fn unit_test() {
    let u0: Unit = Default::default();
    assert_eq!(u0.clone().to_string(), "0");

    let u1: Unit = Unit {
        coef: BigInt::from(3),
        .. Default::default()
    };
    assert_eq!(u1.clone().to_string(), "3");
}


