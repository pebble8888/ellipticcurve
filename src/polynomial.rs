extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;
//use num_traits::One;
use num_traits::Zero;
//use std::cmp::Ordering;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Sub;
use std::ops::SubAssign;
use std::ops::Mul;
//use std::ops::MulAssign;
use std::ops::Div;
use std::ops::Rem;
use std::ops::Neg;
use std::fmt;
use super::unit;

type Unit = unit::Unit;
//type UnitBuilder = unitbuilder::UnitBuilder;

#[derive(Debug, Clone)]
pub struct Polynomial {
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

impl AddAssign for Polynomial {
    fn add_assign(&mut self, other: Self) {
        let mut other_units = other.units.clone();
        self.units.append(&mut other_units);
        self.normalize();
    }
}

impl Sub for Polynomial {
    type Output = Polynomial;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl SubAssign for Polynomial {
    fn sub_assign(&mut self, other: Self) {
        let minus_other = - other;
        let mut tmp_units = minus_other.units.clone();
        self.units.append(&mut tmp_units);
        self.normalize();
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
        if other.is_zero() {
            panic!("other.is_zero()");
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

impl Rem for Polynomial {
    type Output = Polynomial;
    fn rem(self, other: Self) -> Self {
        assert!(!self.has_y());
        assert!(!other.has_y());
        let o_h = other.highest_unit_x();
        let mut tmp = self.clone();
        loop {
            let tmp_h = tmp.highest_unit_x();
            if tmp_h.xpow < o_h.xpow {
                break;
            }
            let q = Polynomial{ units: vec![
                Unit { 
                    coef: tmp_h.coef.clone() / o_h.coef.clone(),
                    xpow: tmp_h.xpow.clone() - o_h.xpow.clone(), 
                    ypow: Zero::zero(),
                },
            ]};
            tmp -= q * other.clone();
        }
        tmp
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = self.clone();
        s.normalize();
        if s.is_zero() {
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
    pub fn new() -> Self {
        Polynomial { units: vec![] }
    }

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
    pub fn power(&self, val: BigInt) -> Self {
        if val == BigInt::from(0) {
            assert!(false);
        }
        let mut m = self.clone();
        for _i in num_iter::range(BigInt::from(0), val-1) {
            m = m.clone() * self.clone();
        }
        m.normalize();
        m
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

    pub fn has_y(&self) -> bool {
        for i in &self.units {
            if i.has_y() {
                return true;
            }
        }
        false
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

    pub fn ec_reduction(&self, a: BigInt, b: BigInt) -> Self {
        let mut t = Polynomial::new();
        for u in &self.units {
            if u.ypow >= BigInt::from(2) {
                let yy = u.ypow.clone() / BigInt::from(2);
                let mut e = Polynomial {
                    units: [u.clone()].to_vec(),
                };
                e = e.clone() / Polynomial { units: [
                Unit {
                    coef: BigInt::from(1), 
                    xpow: BigInt::from(0),
                    ypow: yy.clone() * 2,
                }].to_vec()};
                e = e.clone() * (Polynomial { units: [
                Unit {
                    coef: BigInt::from(1),
                    xpow: BigInt::from(3),
                    ypow: BigInt::from(0),
                }, 
                Unit {
                    coef: a.clone(),
                    xpow: BigInt::from(1),
                    ypow: BigInt::from(0),
                },
                Unit {
                    coef: b.clone(),
                    xpow: BigInt::from(0),
                    ypow: BigInt::from(0),
                },
                ].to_vec() }.power(yy.clone()));
                t += e;
            } else {
                let e = Polynomial {
                    units: [u.clone()].to_vec(),
                };
                t += e;
            }
        }
        t.normalize();
        t
    }

    //pub fn remainder(&self, pol: Polynomial, prime: BigInt) -> Self {
    //
    //}

    // TODO:
    //pub fn expmod(&self, val: usize) -> Self {
    //}
}

#[test]
fn polynmomial_test() { 

    use num_traits::One;
    use super::unitbuilder;
    type UnitBuilder = unitbuilder::UnitBuilder;

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

    // power
    let p39 = p38.clone().power(BigInt::from(2));
    assert_eq!(p39.to_string(), "4 x^4 y^8");

    let p40 = p37.clone().power(BigInt::from(2));
    assert_eq!(p40.to_string(), "25 x^6 y^10 + 10 x^5 y^9 + x^4 y^8");

    // ec_reduction
    let u41 = Unit {
        coef: BigInt::from(1),
        xpow: BigInt::from(0),
        ypow: BigInt::from(2),
    };
    let p41 = Polynomial {
        units: vec![u41],
    };
    let p42 = p41.clone().ec_reduction(BigInt::from(2), BigInt::from(7));
    assert_eq!(p42.to_string(), "x^3 + 2 x + 7");

    assert_eq!(p42.highest_unit_x().to_string(), "x^3");

    // rem
    let u43 = UnitBuilder::new().coef(1).xpow(4).finalize();
    let u44 = UnitBuilder::new().coef(2).xpow(1).finalize();
    let p47 = Polynomial { units: vec![u43, u44] };

    let u45 = UnitBuilder::new().coef(1).xpow(2).finalize();
    let u46 = UnitBuilder::new().coef(3).finalize();
    let p48 = Polynomial { units: vec![u45, u46] };
    assert_eq!(p47.clone().to_string(), "x^4 + 2 x");
    assert_eq!(p48.clone().to_string(), "x^2 + 3");

    assert_eq!((p47 % p48).to_string(), "2 x + 9"); 
}

