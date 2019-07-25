extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use num_traits::pow;
use std::cmp::Ordering;
use std::ops::Add;
use std::ops::Sub;
use std::ops::Mul;
use std::ops::Div;
//use std::ops::Rem;
use std::ops::Neg;
use std::fmt;
use super::unit;

type Unit = unit::Unit;

#[derive(Debug, Clone)]
struct Polynomial {
    pub units: Vec<Unit>,
}

impl Add for Polynomial {
    type Output = Polynomial;
    fn add(self, other: Self) -> Self {
        let mut units = self.units;
        let mut other_units = other.units;
        units.append(&mut other_units);
        let mut pol = Self { units: units };
        pol.normalize();
        pol
    }
}

impl Sub for Polynomial {
    type Output = Polynomial;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl Neg for Polynomial {
    type Output = Polynomial;
    fn neg(self) -> Polynomial {
        let mut units: Vec<Unit> = Vec::new();
        for u in &self.units {
            units.push(-u.clone());
        }
        Polynomial { units: units }
    }
}

impl Mul for Polynomial {
    type Output = Polynomial;
    fn mul(self, other: Self) -> Self {
        let mut units: Vec<Unit> = Vec::new();
        for i in &self.units {
            for j in &other.units {
                let u: Unit = i.clone() * j.clone();
                units.push(u);
            }
        }
        let mut pol = Polynomial {
            units: units,
        };
        pol.normalize();
        pol
    }
}

impl Div for Polynomial {
    type Output = Polynomial;
    fn div(self, other: Self) -> Self {
        if other.units.len() == 0 {
            panic!("other.units.len() == 0");
        }
        else if other.units.len() == 1 {
            let mut units: Vec<Unit> = Vec::new();
            for i in &self.units {
                let u = i.clone() / other.units[0].clone();
                units.push(u);
            }
            let mut pol = Polynomial {
                units: units,
            };
            pol.normalize();
            pol
        } else {
            panic!("other.units.len() >= 2");
        }
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = self.clone();
        s.normalize();
        if s.units.len() == 0 {
            return write!(f, "0");
        }
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
        while i < units.len() {
            if units[i].coef.clone() == BigInt::from(0) {
                units.remove(i);
            } else if i+1 < units.len() && units[i].equal_order(&units[i+1]) {
                let coef = units[i].coef.clone() + units[i+1].coef.clone();
                if coef.clone() != BigInt::from(0) {
                    let mut unit = units[i].clone();
                    unit.coef = coef.clone();
                    units.insert(i, unit);
                    units.remove(i+1);
                    units.remove(i+1);
                    i = i + 1;
                } else {
                    units.remove(i);
                    units.remove(i);
                }
            } else {
                i = i + 1;
            }
        }
    }

    pub fn modular(&self, val: BigInt) -> Self {
        let mut units: Vec<Unit> = Vec::new();
        for i in &self.units {
            let u = i.modular(val.clone());
            units.push(u);
        }
        let mut pol = Polynomial {
            units: units,
        };
        pol.normalize();
        pol
    }

    pub fn highest_unit_x(&self) -> Unit {
        let mut pol = self.clone();
        pol.normalize(); 
        pol.units[0].clone()
    }

    pub fn is_zero(&self) -> bool {
        let mut pol = self.clone();
        pol.normalize();
        pol.units.len() == 0
    }

    pub fn to_frob(&self, n: usize) -> Self {
        let mut units: Vec<Unit> = Vec::new();
        for u in &self.units {
            units.push(u.to_frob(n));
        }
        let mut pol = Polynomial {
            units: units,
        };
        pol.normalize();
        pol
    }

    /*
    pub fn ec_reduction(&self, a: BigInt, b: BigInt) -> Self {
        let mut t = Polynomial {
            units: [].to_vec(),
        };
        for u in &self.units {
            if u.ypow >= BigInt::from(2) {
                let yy = u.ypow.modular(BigInt::from(2));
                let e = Polynomial {
                }
            }
        }
    }
    */

    //pub fn remainder(&self, pol: Polynomial, prime: BigInt) -> Self {
    //
    //}

    // TODO:
    //pub fn expmod(&self, val: usize) -> Self {
    //}
}

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
    let p4 = p2.clone() + p3;
    assert_eq!(p4.to_string(), "x^4 y^2 + 3 x^2 y^4");

    let u6 = u1.clone() * u1.clone();
    assert_eq!(u6.to_string(), "x^8 y^4");
    
    let p5 = p1.clone() * p4.clone();
    assert_eq!(p5.to_string(), "x^8 y^4 + 6 x^6 y^6 + 9 x^4 y^8");

    // Neg
    let u7 = u1.clone();
    let u8 = - u7; 
    assert_eq!(u8.to_string(), "- x^4 y^2"); 

    let p8 = p1.clone();
    let p9 = - p8; 
    assert_eq!(p9.to_string(), "- x^4 y^2 - 3 x^2 y^4");

    // Sub
    let p10 = p1.clone() - p2.clone(); 
    assert_eq!(p10.to_string(), "3 x^2 y^4");

    let p11 = p5.clone() - p5.clone();
    assert_eq!(p11.to_string(), "0");

    let u21 = u2.clone().power(3);
    assert_eq!(u21.to_string(), "27 x^6 y^12");

    // Modular
    let u33 = Unit {
        coef: BigInt::from(37),
        xpow: BigInt::from(2),
        ypow: BigInt::from(4),
    };

    let u35 = Unit {
        coef: BigInt::from(35),
        xpow: BigInt::from(3),
        ypow: BigInt::from(5),
    };

    let u34 = u33.clone().modular(BigInt::from(24));
    assert_eq!(u34.to_string(), "13 x^2 y^4");

    let p36 = Polynomial {
        units: vec![u33.clone(), u35.clone()],
    };
    let p37 = p36.clone().modular(BigInt::from(6));
    assert_eq!(p37.to_string(), "5 x^3 y^5 + x^2 y^4");

    let p38 = p36.clone().modular(BigInt::from(5));
    assert_eq!(p38.to_string(), "2 x^2 y^4");

}

