extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;
use num_traits::Zero;
use num_traits::One;
use std::fmt;
use std::ops;
use super::unit;
use super::bigint::{Inverse, DivFloor};

type Unit = unit::Unit;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub units: Vec<Unit>,
}

impl_op_ex!(+ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut units = a.units.clone();
    let mut other_units = b.units.clone();
    units.append(&mut other_units);
    let mut pol = Polynomial { units: units };
    pol.normalize();
    pol
});

impl_op_ex!(+= |a: &mut Polynomial, b: &Polynomial| {
    a.units.append(&mut b.units.clone());
    a.normalize();
});

impl_op_ex!(- |a: &Polynomial, b: &Polynomial| -> Polynomial {
    a.clone() + (-b.clone())
});

impl_op_ex!(-= |a: &mut Polynomial, b: &Polynomial| {
    a.units.append(&mut (-b).units.clone());
    a.normalize();
});

impl_op_ex!(- |a: &Polynomial| -> Polynomial {
    let mut units: Vec<Unit> = Vec::with_capacity(a.units.len());
    for u in &a.units {
        units.push(-u.clone());
    }
    Polynomial { units: units }
});

impl_op_ex!(* |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut units: Vec<Unit> = Vec::with_capacity(a.units.len() * b.units.len());
    for i in &a.units {
        for j in &b.units {
            let u: Unit = i * j;
            units.push(u);
        }
    }
    let mut pol = Polynomial {
        units: units,
    };
    pol.normalize();
    pol
});

impl_op_ex!(*= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() * b;
    a.units.clear();
    a.units = c.units; 
    a.normalize();
});

impl_op_ex!(/ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    if b.is_zero() {
        panic!("b.is_zero()");
    }
    else if b.units.len() == 1 {
        let mut units: Vec<Unit> = Vec::with_capacity(a.units.len());
        for i in &a.units {
            let u = i / &b.units[0];
            units.push(u);
        }
        let mut pol = Polynomial {
            units: units,
        };
        pol.normalize();
        pol
    } else {
        panic!("b.units.len() >= 2");
    }
});

impl_op_ex!(/= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() / b;
    a.units.clear();
    a.units = c.units;
    a.normalize();
});

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
        Polynomial { units: vec![Unit::new()] }
    }

    pub fn normalize(&mut self) {
        println!("{} normalize()", line!());
        let units = &mut self.units;
        units.sort_by(|a, b| b.cmp(a));
        let mut i = 0; 
        while i < units.len() {
            if units[i].coef == Zero::zero() {
                units.remove(i);
            } else if i+1 < units.len() && units[i].equal_order(&units[i+1]) {
                let coef = &units[i].coef + &units[i+1].coef;
                if coef != Zero::zero() {
                    let mut unit = units[i].clone();
                    unit.coef = coef;
                    units.insert(i, unit);
                    units.remove(i+1);
                    units.remove(i+1);
                } else {
                    units.remove(i);
                    units.remove(i);
                }
            } else {
                i += 1;
            }
        }
        println!("{} normalize()", line!());
    }

    pub fn power(&self, n: &BigInt) -> Self {
        println!("{} {}", line!(), n);
        assert!(*n >= Zero::zero(), "n:{}", n.to_string());
        if *n == Zero::zero() {
            return Unit {
                coef: BigInt::from(1),
                xpow: Zero::zero(),
                ypow: Zero::zero(),
            }.to_pol();
        }
        
        let mut m = self.clone();
        for _i in num_iter::range(BigInt::from(0), n-1) {
            m *= self;
        }
        m.normalize();
        println!("{}", line!());
        m
    }

    pub fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self {
        println!("{}: n:{}", line!(), n);
        assert!(*n >= Zero::zero());
        let mut b = self.clone();
        b = b.modular(p);
        let mut r = Unit {
                coef: BigInt::from(1),
                xpow: Zero::zero(),
                ypow: Zero::zero(),
        }.to_pol();
        let mut e = n.clone();
        while &e > &One::one() {
            println!("{} e:{}", line!(), e);
            if &e % 2 != Zero::zero() {
                r *= &b;
                r = r.modular(p);
            }
            b *= b.clone();
            b = b.modular(p);
            e /= 2;
        }
        //println!("{} r:{} b:{}", line!(), r.to_string(), b.to_string());
        r *= b;
        //println!("{}", line!());
        r = r.modular(p);
        //println!("{}", line!());
        r
    }

    pub fn polynomial_modular(&self, other: &Polynomial, p: &BigInt) -> Self {
        assert!(!self.has_y(), "{}", self.to_string());
        assert!(!other.has_y(), "{}", other.to_string());
        let h = other.highest_unit_x();
        let mut tmp = self.clone();
        loop {
            let tmp_h = tmp.highest_unit_x();
            if tmp_h.xpow < h.xpow {
                break;
            }
            let q = Unit { 
                        coef: &tmp_h.coef * &h.coef.inverse(p),
                        xpow: &tmp_h.xpow - &h.xpow, 
                        ypow: Zero::zero(),
                    }.to_pol();
            tmp -= q * other;
            tmp = tmp.modular(&p);
        }
        tmp.modular(&p)
    }

    pub fn modular(&self, p: &BigInt) -> Self {
        let mut tmp = self.clone();
        tmp.normalize();
        let mut units: Vec<Unit> = Vec::with_capacity(tmp.units.len());
        for j in &tmp.units {
            let u = j.modular(p);
            if p != &Zero::zero() {
                assert!(&u.coef < p, "u.coef {}", &u.coef);
            }
            units.push(u);
        }
        let mut pol = Polynomial {
            units: units,
        };
        pol.normalize();

        for i in &pol.units {
            if p != &Zero::zero() {
                assert!(&i.coef < p, "i.coef {} {}", &i.coef, pol.to_string());
            }
        }
        pol
    }

    pub fn highest_unit_x(&self) -> Unit {
        let mut pol = self.clone();
        pol.normalize(); 
        if pol.is_zero() {
            return Unit::new();
        }
        pol.units[0].clone()
    }

    pub fn is_zero(&self) -> bool {
        let mut pol = self.clone();
        pol.normalize();
        pol.units.len() == 0
    }

    pub fn is_gcd_one(&self, other: &Self, p: &BigInt) -> bool {
        let s = self.highest_unit_x();
        let o = other.highest_unit_x();
        if s.xpow < o.xpow {
            let m = other.polynomial_modular(self, p);
            return !m.is_zero();
        }
        let m = self.polynomial_modular(other, p);
        return !m.is_zero();
    }

    pub fn has_y(&self) -> bool {
        for i in &self.units {
            if i.has_y() {
                return true;
            }
        }
        false
    }

    pub fn to_frob(&self, n: &BigInt) -> Self {
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

    pub fn to_y_power(&self, n: &BigInt) -> Self {
        let mut units: Vec<Unit> = Vec::new();
        for u in &self.units {
            units.push(u.to_y_power(n));
        }
        let mut pol = Polynomial {
            units: units,
        };
        pol.normalize();
        pol
    }

    pub fn ec_reduction(&self, a: &BigInt, b: &BigInt) -> Self {
        let mut t = Polynomial::new();
        for u in &self.units {
            if u.ypow >= BigInt::from(2) {
                let yy = u.ypow.clone().div_floor(&BigInt::from(2));
                let mut e = u.clone().to_pol();
                e /= Polynomial { units: [
                Unit {
                    coef: BigInt::from(1), 
                    xpow: BigInt::from(0),
                    ypow: yy.clone() * 2,
                }].to_vec()};
                let ee = Polynomial { units: [
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
                ].to_vec() };
                // power() is faster than power_modular()
                println!("{}", line!());
                let ee = ee.power(&yy.clone());
                println!("{}", line!());
                e *= ee;
                println!("{}", line!());
                t += e;
            } else {
                //println!("{}", line!());
                t += u.clone().to_pol();
            }
        }
        println!("{}", line!());
        t.normalize();
        println!("{}", line!());
        t
    }
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
    let p4 = &p2 + &p3;
    assert_eq!(p4.to_string(), "x^4 y^2 + 3 x^2 y^4");

    let u6 = &u1 * &u1;
    assert_eq!(u6.to_string(), "x^8 y^4");
    
    let p5 = &p1 * &p4;
    assert_eq!(p5.to_string(), "x^8 y^4 + 6 x^6 y^6 + 9 x^4 y^8");

    // Neg
    let u8 = - &u1;
    assert_eq!(u8.to_string(), "- x^4 y^2"); 

    let p9 = - &p1; 
    assert_eq!(p9.to_string(), "- x^4 y^2 - 3 x^2 y^4");

    // Sub
    let p10 = &p1 - &p2; 
    assert_eq!(p10.to_string(), "3 x^2 y^4");

    let p11 = &p5 - &p5;
    assert_eq!(p11.to_string(), "0");

    let u21 = u2.power(&BigInt::from(3));
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

    let u34 = u33.modular(&BigInt::from(24));
    assert_eq!(u34.to_string(), "13 x^2 y^4");

    let p36 = Polynomial {
        units: vec![u33.clone(), u35.clone()],
    };
    let p37 = p36.modular(&BigInt::from(6));
    assert_eq!(p37.to_string(), "5 x^3 y^5 + x^2 y^4");

    let p38 = p36.modular(&BigInt::from(5));
    assert_eq!(p38.to_string(), "2 x^2 y^4");

    // power
    let p39 = p38.power(&BigInt::from(2));
    assert_eq!(p39.to_string(), "4 x^4 y^8");

    let p3_39 = p38.power(&BigInt::from(3));
    assert_eq!(p3_39.to_string(), "8 x^6 y^12");

    let p40 = p37.power(&BigInt::from(2));
    assert_eq!(p40.to_string(), "25 x^6 y^10 + 10 x^5 y^9 + x^4 y^8");

    // ec_reduction
    let u41 = Unit {
        coef: BigInt::from(1),
        xpow: BigInt::from(1),
        ypow: BigInt::from(3),
    };
    let p41 = Polynomial {
        units: vec![u41],
    };
    let p42 = p41.ec_reduction(&BigInt::from(2), &BigInt::from(7));
    assert_eq!(p42.to_string(), "x^4 y + 2 x^2 y + 7 x y");

    assert_eq!(p42.highest_unit_x().to_string(), "x^4 y");

    // rem
    let u43 = UnitBuilder::new().coef(3).xpow(4).finalize();
    let u44 = UnitBuilder::new().coef(2).xpow(1).finalize();
    let p47 = Polynomial { units: vec![u43, u44] };

    let u45 = UnitBuilder::new().coef(2).xpow(2).finalize();
    let u46 = UnitBuilder::new().coef(3).finalize();
    let p48 = Polynomial { units: vec![u45, u46] };
    assert_eq!(p47.to_string(), "3 x^4 + 2 x");
    assert_eq!(p48.to_string(), "2 x^2 + 3");

    // plynomial_modular
    assert_eq!(p47.polynomial_modular(&p48, &BigInt::from(19)).to_string(), "2 x + 2"); 

    // *=
    let mut p49 = p47.clone();
    p49 *= &p48;
    assert_eq!(p49.to_string(), "6 x^6 + 9 x^4 + 4 x^3 + 6 x");

    // /=
    let u50 = UnitBuilder::new().coef(6).xpow(6).finalize();
    let u51 = UnitBuilder::new().coef(12).xpow(4).finalize();
    let u52 = UnitBuilder::new().coef(3).xpow(2).finalize();
    let mut p50 = &u50 + &u51;
    let p51 = u52.to_pol();
    p50 /= &p51;
    assert_eq!(p50.to_string(), "2 x^4 + 4 x^2");

}

