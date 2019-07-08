extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use std::cmp::Ordering;
use std::ops::Add;

#[derive(Debug, Clone)]
struct Polynomial {
    units: Vec<Unit>,
}

// 3 x^4 y^5
#[derive(Debug, Clone, PartialEq, Eq)]
struct Unit {
    coef: BigInt, // 3
    xpow: BigInt, // 4
    ypow: BigInt, // 5
}

use std::fmt;

// Unit

impl Add for Polynomial {
    type Output = Polynomial;
    fn add(self, other: Self) -> Self {
        let mut units = self.units;
        let mut other_units = other.units;
        units.append(&mut other_units);
        let mut s = Self { units: units };
        s.normalize();
        s
    }
}

impl Unit {
    pub fn equal_order(&self, other: &Self) -> bool {
        return &self.xpow == &other.xpow && self.ypow == other.ypow 
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

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = self.clone();
        s.normalize();
        let mut st = String::new();
        let mut i = 0;
        for u in &s.units {
            if u.coef > Zero::zero() && i != 0 {
                st.push_str("+ ");
            }
            st.push_str(&u.to_string());
            st.push_str(" ");
            i = i + 1;
        }
        write!(f, "{}", st.trim_end())
    }
}

impl Polynomial {
    pub fn normalize(&mut self) {
        let units = &mut self.units;
        units.sort_by(|a, b| b.cmp(a));
        let mut i = 0; 
        while i+1 < units.len() {
            if units[i].equal_order(&units[i+1]) {
                let coef = units[i].coef.clone() + units[i+1].coef.clone();
                let mut unit = units[i].clone();
                unit.coef = coef;
                units.insert(i, unit);
                units.remove(i+1);
                units.remove(i+1);
            }
            i = i + 1;
        }
    }
}

impl fmt::Display for Unit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.coef == One::one() {
            if self.xpow == One::one() && self.ypow == One::one() {
                write!(f, "1")
            } else {
                let mut st = String::new();
                if self.xpow != One::one() {
                    st.push_str("x^");
                    st.push_str(&self.xpow.to_string());
                    st.push_str(" ");
                }
                if self.ypow != One::one() {
                    st.push_str("y^");
                    st.push_str(&self.ypow.to_string());
                }
                write!(f, "{}", st.trim_end())
            }
        } else {
            let mut st = String::new();
            st.push_str(&self.coef.to_string());
            st.push_str(" ");
            if self.xpow != One::one() {
                st.push_str("x^");
                st.push_str(&self.xpow.to_string());
                st.push_str(" ");
            }
            if self.ypow != One::one() {
                st.push_str("y^");
                st.push_str(&self.ypow.to_string());
            }
            write!(f, "{}", st.trim_end())
        }
    }
}

// test unit

#[test]
fn unit() {
    let u1 = Unit {
        coef: One::one(),
        xpow: BigInt::from(4),
        ypow: BigInt::from(2),
    };
    let u2 = Unit {
        coef: BigInt::from(3),
        xpow: BigInt::from(2),
        ypow: BigInt::from(4),
    };
    assert_eq!(u1.to_string(), "x^4 y^2");
    assert_eq!(u2.to_string(), "3 x^2 y^4");

    let p1 = Polynomial {
        units: vec![u1.clone(), u2.clone()],
    };
    assert_eq!(p1.to_string(), "x^4 y^2 + 3 x^2 y^4");

    let p2 = Polynomial {
        units: vec![u1.clone()],
    };
    let p3 = Polynomial {
        units: vec![u2.clone()],
    };
    let p4 = p2 + p3;
    assert_eq!(p4.to_string(), "x^4 y^2 + 3 x^2 y^4");
}

