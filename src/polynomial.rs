use num_integer::Integer;
use std::collections::BTreeMap;
use num_bigint::BigInt;
use num_traits::Zero;
use num_traits::One;
use std::fmt;
use std::ops;
use super::unit;
use super::bigint::{Inverse};
use super::unitbuilder;

type Unit = unit::Unit;
type UnitBuilder = unitbuilder::UnitBuilder;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub units: BTreeMap<unit::UnitKey, BigInt>
}

// +
impl_op_ex!(+ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut pol = a.clone();
    for (bk, bv) in &b.units {
        if let Some(av) = pol.units.get_mut(&bk) {
            *av += bv;
            if av.is_zero() {
                pol.units.remove(&bk.clone());
            }
        } else {
            pol.units.insert(bk.clone(), bv.clone());
        }
    }
    pol
});

impl_op_ex!(+ |a: &Polynomial, b: &Unit| -> Polynomial {
    a + b.to_pol()
});

impl_op_ex!(+ |a: &Unit, b: &Polynomial| -> Polynomial {
    a.to_pol() + b
});

// +=
impl_op_ex!(+= |a: &mut Polynomial, b: &Polynomial| {
    for (bk, bv) in &b.units {
        if let Some(av) = a.units.get_mut(&bk) {
            *av += bv;
            if av.is_zero() {
                a.units.remove(&bk.clone());
            }
        } else {
            a.units.insert(bk.clone(), bv.clone());
        }
    }
});

impl_op_ex!(+= |a: &mut Polynomial, b: &Unit| {
    if let Some(av) = a.units.get_mut(&b.key) {
        *av += &b.coef;
        if av.is_zero() {
            a.units.remove(&b.key.clone());
        }
    } else {
        a.units.insert(b.key.clone(), b.coef.clone());
    }
});

// -
impl_op_ex!(- |a: &Polynomial, b: &Polynomial| -> Polynomial {
    a.clone() + (-b.clone())
});

impl_op_ex!(- |a: &Polynomial, b: &Unit | -> Polynomial {
    a.clone() + (-b.to_pol())
});

impl_op_ex!(- |a: &Unit, b: &Polynomial| -> Polynomial {
    a.to_pol() + (-b.clone())
});

// -=
impl_op_ex!(-= |a: &mut Polynomial, b: &Polynomial| {
    for (bk, bv) in &b.units {
        if let Some(av) = a.units.get_mut(&bk) {
            *av -= bv;
            if av.is_zero() {
                a.units.remove(&bk.clone());
            }
        } else {
            a.units.insert(bk.clone(), - bv.clone());
        }
    }
});

// Neg
impl_op_ex!(- |a: &Polynomial| -> Polynomial {
    let mut pol = Polynomial::new();
    for (k, coef) in &a.units {
        pol.units.insert(k.clone(), -coef);
    }
    pol
});

// *
impl_op_ex!(* |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut pol = Polynomial::new();
    for (ik, iv) in &a.units {
        let i = Unit::from(ik, iv);
        for (jk, jv) in &b.units {
            let j = Unit::from(jk, jv);
            let l = i.clone() * j;
            if let Some(lv) = pol.units.get_mut(&l.key) {
                *lv += &l.coef;
            } else {
                pol.units.insert(l.key, l.coef);
            }
        }
    }
    pol
});

impl_op_ex!(* |a: &Polynomial, b: &Unit| -> Polynomial {
    a * b.to_pol()
});

impl_op_ex!(* |a: &Unit, b: &Polynomial| -> Polynomial {
    a.to_pol() * b
});

// *=
impl_op_ex!(*= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() * b;
    a.units.clear();
    a.units = c.units; 
});

// /
impl_op_ex!(/ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    if b.is_zero() {
        panic!("b.is_zero()");
    }
    else if b.units.len() == 1 {
        let mut pol = Polynomial::new();
        for (ak, av) in &a.units {
            let i = Unit::from(ak, av);
            let (k, v) = b.units.iter().next().unwrap();
            let u2 = Unit::from(k, v); 
            let u = i / &u2;
            pol.units.insert(u.key, u.coef);
        }
        pol
    } else {
        panic!("b.units.len() >= 2");
    }
});

impl_op_ex!(/ |a: &Polynomial, b: &Unit| -> Polynomial {
    a / b.to_pol()
});

impl_op_ex!(/ |a: &Unit, b: &Polynomial| -> Polynomial {
    a.to_pol() / b
});

// /=
impl_op_ex!(/= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() / b;
    a.units.clear();
    a.units = c.units;
});

impl_op_ex!(/= |a: &mut Polynomial, b: &Unit| {
    let c = a.clone() / b.to_pol();
    a.units.clear();
    a.units = c.units;
});

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let  s = self.clone();
        if s.is_zero() {
            return write!(f, "0");
        }
        let mut st = String::new();
        let mut i = 0;
        for (k, coef) in s.units.iter().rev() {
            if coef > &BigInt::from(0) && i != 0 {
                st.push_str("+ ");
            }
            st.push_str(&Unit::from(k, coef).to_string());
            st.push_str(" ");
            i = i + 1;
        }
        write!(f, "{}", st.trim_end())
    }
}

impl Polynomial {
    pub fn new() -> Self {
        Polynomial {
            units: BTreeMap::new(),  
        }
    }

    pub fn power(&self, n: &BigInt) -> Self {
        assert!(n >= &Zero::zero(), "n:{}", n.to_string());
        if n == &Zero::zero() {
            UnitBuilder::new().coef_i(1).finalize().to_pol();
        }
        
        let mut m = self.clone();
        for _i in num_iter::range(BigInt::from(0), n-1) {
            m *= self;
        }
        m
    }

    pub fn power_i(&self, n: i64) -> Self {
        let n = BigInt::from(n);
        self.power(&n)
    }

    pub fn square(&self) -> Self {
        self.power(&BigInt::from(2))
    }

    pub fn cube(&self) -> Self {
        self.power(&BigInt::from(3))
    }

    pub fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self {
        assert!(*n >= Zero::zero());
        let mut b = self.modular(p);
        let mut r = UnitBuilder::new().coef_i(1).finalize().to_pol();
        let mut e = n.clone();
        while &e > &One::one() {
            if &e % 2 != Zero::zero() {
                r *= &b;
                r = r.modular(p);
            }
            b *= b.clone();
            b = b.modular(p);
            e /= 2;
        }
        r *= b;
        r = r.modular(p);
        r
    }

    pub fn polynomial_modular(&self, other: &Polynomial, p: &BigInt) -> Self {
        assert!(!self.has_y(), "!self.has_y()");
        assert!(!other.has_y(), "!other.has_y()");
        let oh = other.highest_unit_x();
        let mut r = self.clone();
        loop {
            let rh = r.highest_unit_x();
            if rh.xpow() < oh.xpow() {
                break;
            }
            let q = UnitBuilder::new()
                    .coef(&(rh.coef.clone() * oh.coef.inverse(p)))
                    .xpow(&(rh.xpow() - oh.xpow()))
                    .finalize().to_pol();
            let q = q.modular(&p);
            let d = (q * other).modular(&p);
            r -= d;
            r = r.modular(&p);
        }
        r.modular(&p)
    }

    pub fn modular(&self, p: &BigInt) -> Self {
        let tmp = self.clone();
        let mut pol = Polynomial::new();  
        for (k, coef) in tmp.units {
            let c = coef.mod_floor(p);
            if p != &Zero::zero() {
                assert!(&c < p, "c {}", &c);
            }
            if !c.is_zero() {
                pol.units.insert(k, c);
            }
        }

        for (_, coef) in &pol.units {
            if p != &Zero::zero() {
                assert!(coef < p, "{} {}", &coef, pol.to_string());
            }
        }
        pol
    }

    pub fn highest_unit_x(&self) -> Unit {
        if self.is_zero() {
            return Unit::new();
        }
        for (k, v) in self.units.iter().rev() {
            return Unit::from(k, v);
        }
        panic!("highest_unit_x assert!");
    }

    pub fn is_zero(&self) -> bool {
        self.units.len() == 0
    }

    pub fn is_gcd_one(&self, other: &Self, p: &BigInt) -> bool {
        let s = self.highest_unit_x();
        let o = other.highest_unit_x();
        if s.xpow() < o.xpow() {
            let m = other.polynomial_modular(self, p);
            return !m.is_zero();
        }
        let m = self.polynomial_modular(other, p);
        return !m.is_zero();
    }

    pub fn has_y(&self) -> bool {
        for (k, v) in &self.units {
            let i = Unit::from(k, v);
            if i.has_y() {
                return true;
            }
        }
        false
    }

    pub fn to_frob(&self, n: &BigInt) -> Self {
        let mut pol = Polynomial::new();
        for (k, v) in &self.units {
            let u = Unit::from(k, v);
            let u = u.to_frob(n);
            pol.units.insert(u.key, u.coef);
        }
        pol
    }

    pub fn to_y_power(&self, n: &BigInt) -> Self {
        let mut pol = Polynomial::new();
        for (k, v) in &self.units {
            let u = Unit::from(k, v);
            let u = u.to_y_power(n);
            pol.units.insert(u.key, u.coef);
        }
        pol
    }

    pub fn reduction(&self, a: &BigInt, b: &BigInt) -> Self {
        let mut t = Polynomial::new();
        for (k, v) in &self.units {
            let u = Unit::from(k, v);
            if u.ypow() >= &BigInt::from(2) {
                let yy = u.ypow().clone().div_floor(&BigInt::from(2));
                let mut e = u.clone().to_pol();
                e /= UnitBuilder::new().coef_i(1).ypow(&(yy.clone() * 2)).finalize();

                let ee = UnitBuilder::new().coef_i(1).xpow_i(3).finalize()
                       + UnitBuilder::new().coef(&a.clone()).xpow_i(1).finalize()
                       + UnitBuilder::new().coef(&b.clone()).finalize();
                // power() is faster than power_modular()
                let ee = ee.power(&yy.clone());
                e *= ee;
                t += e;
            } else {
                t += u.clone();
            }
        }
        t
    }

    pub fn reduction_modular(&self, a: &BigInt, b: &BigInt, p: &BigInt) -> Self {
        let mut t = Polynomial::new();
        for (k, v) in &self.units {
            let u = Unit::from(k, v);
            if u.ypow() >= &BigInt::from(2) {
                let yy = u.ypow().clone().div_floor(&BigInt::from(2));
                let mut e = u.clone().to_pol();
                e /= UnitBuilder::new().coef_i(1).ypow(&(yy.clone() * 2)).finalize();

                let ee = UnitBuilder::new().coef_i(1).xpow_i(3).finalize()
                       + UnitBuilder::new().coef(&a.clone()).xpow_i(1).finalize()
                       + UnitBuilder::new().coef(&b.clone()).finalize();
                // power() is faster than power_modular()
                let ee = ee.power_modular(&yy.clone(), p);
                let ee = ee.reduction_modular(a, b, p);
                e *= ee;
                t += e;
            } else {
                t += u.clone();
            }
            t = t.modular(p);
        }
        t
    }
}

#[test]
fn polynmomial_test() { 

    use super::unitbuilder;
    type UnitBuilder = unitbuilder::UnitBuilder;

    let u1 = UnitBuilder::new().coef_i(1).xpow_i(4).ypow_i(2).finalize();
    let u2 = UnitBuilder::new().coef_i(3).xpow_i(2).ypow_i(4).finalize();
    assert_eq_str!(u1, "x^4 y^2");
    assert_eq_str!(u2, "3 x^2 y^4");

    let p1 = u1.clone() + u2.clone();
    assert_eq_str!(p1, "x^4 y^2 + 3 x^2 y^4");

    let p2 = u1.to_pol();
    let p3 = u2.to_pol();
    let p4 = &p2 + &p3;
    assert_eq_str!(p4, "x^4 y^2 + 3 x^2 y^4");

    let u6 = &u1 * &u1;
    assert_eq_str!(u6, "x^8 y^4");
    
    let p5 = &p1 * &p4;
    assert_eq_str!(p5, "x^8 y^4 + 6 x^6 y^6 + 9 x^4 y^8");

    // Neg
    let p9 = - &p1; 
    assert_eq_str!(p9, "- x^4 y^2 - 3 x^2 y^4");

    // Sub
    let p10 = &p1 - &p2; 
    assert_eq_str!(p10, "3 x^2 y^4");

    let p11 = &p5 - &p5;
    assert_eq_str!(p11, "0");

    let u21 = u2.power_i(3);
    assert_eq_str!(u21, "27 x^6 y^12");

    // Modular
    let u33 = UnitBuilder::new().coef_i(37).xpow_i(2).ypow_i(4).finalize();

    let u35 = UnitBuilder::new().coef_i(35).xpow_i(3).ypow_i(5).finalize();

    let u34 = u33.modular(&BigInt::from(24));
    assert_eq_str!(u34, "13 x^2 y^4");

    let p36 = u33 + u35;
    let p37 = p36.modular(&BigInt::from(6));
    assert_eq_str!(p37, "5 x^3 y^5 + x^2 y^4");

    let p38 = p36.modular(&BigInt::from(5));
    assert_eq_str!(p38, "2 x^2 y^4");

    //println!("{}", line!());

    // power
    let p39 = p38.power_i(2);
    assert_eq_str!(p39, "4 x^4 y^8");

    let p3_39 = p38.power_i(3);
    assert_eq_str!(p3_39, "8 x^6 y^12");

    let p40 = p37.power_i(2);
    assert_eq_str!(p40, "25 x^6 y^10 + 10 x^5 y^9 + x^4 y^8");

    // reduction
    let u41 = UnitBuilder::new().coef_i(1).xpow_i(1).ypow_i(3).finalize();
    let p41 = u41.to_pol();
    let p42 = p41.reduction(&BigInt::from(2), &BigInt::from(7));
    assert_eq_str!(p42, "x^4 y + 2 x^2 y + 7 x y");

    let p422 = UnitBuilder::new().coef_i(2).xpow_i(3).finalize().to_pol();
    assert_eq_str!(p422.highest_unit_x(), "2 x^3");

    // rem
    let u43 = UnitBuilder::new().coef_i(3).xpow_i(4).finalize();
    let u44 = UnitBuilder::new().coef_i(2).xpow_i(1).finalize();
    let p47 = u43 + u44;

    let u45 = UnitBuilder::new().coef_i(2).xpow_i(2).finalize();
    let u46 = UnitBuilder::new().coef_i(3).finalize();
    let p48 = u45 + u46;
    assert_eq_str!(p47, "3 x^4 + 2 x");
    assert_eq_str!(p48, "2 x^2 + 3");

    // plynomial_modular
    assert_eq_str!(p47.polynomial_modular(&p48, &BigInt::from(19)), "2 x + 2"); 

    // *=
    let mut p49 = p47.clone();
    p49 *= &p48;
    assert_eq_str!(p49, "6 x^6 + 9 x^4 + 4 x^3 + 6 x");

    // -=
    let mut p80 = UnitBuilder::new().coef_i(6).xpow_i(6).finalize()
                + UnitBuilder::new().coef_i(12).xpow_i(4).finalize();
    let p81 = UnitBuilder::new().coef_i(2).xpow_i(6).finalize().to_pol();
    p80 -= p81;
    assert_eq_str!(p80, "4 x^6 + 12 x^4");

    // /=
    let u50 = UnitBuilder::new().coef_i(6).xpow_i(6).finalize();
    let u51 = UnitBuilder::new().coef_i(12).xpow_i(4).finalize();
    let u52 = UnitBuilder::new().coef_i(3).xpow_i(2).finalize();
    let mut p50 = &u50 + &u51;
    let p51 = u52.to_pol();
    p50 /= &p51;
    assert_eq_str!(p50, "2 x^4 + 4 x^2");
}

