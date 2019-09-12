use num_integer::Integer;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use std::{fmt, ops};
use super::bigint::{Inverse, Power};
use super::term;
use super::termbuilder;
use super::termbuilder::TermBuildable;

type Term = term::Term;
type TermBuilder = termbuilder::TermBuilder;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub terms: BTreeMap<term::Monomial, BigInt>
}

// +
impl_op_ex!(+ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut pol = a.clone();
    for (bk, bv) in &b.terms {
        if let Some(av) = pol.terms.get_mut(&bk) {
            *av += bv;
            if av.is_zero() {
                pol.terms.remove(&bk.clone());
            }
        } else {
            pol.terms.insert(bk.clone(), bv.clone());
        }
    }
    pol
});

impl_op_ex!(+ |a: &Polynomial, b: &Term| -> Polynomial {
    a + b.to_pol()
});

impl_op_ex!(+ |a: &Term, b: &Polynomial| -> Polynomial {
    a.to_pol() + b
});

// +=
impl_op_ex!(+= |a: &mut Polynomial, b: &Polynomial| {
    for (bk, bv) in &b.terms {
        if let Some(av) = a.terms.get_mut(&bk) {
            *av += bv;
            if av.is_zero() {
                a.terms.remove(&bk.clone());
            }
        } else {
            a.terms.insert(bk.clone(), bv.clone());
        }
    }
});

impl_op_ex!(+= |a: &mut Polynomial, b: &Term| {
    if let Some(av) = a.terms.get_mut(&b.monomial) {
        *av += &b.coef;
        if av.is_zero() {
            a.terms.remove(&b.monomial.clone());
        }
    } else {
        a.terms.insert(b.monomial.clone(), b.coef.clone());
    }
});

// -
impl_op_ex!(- |a: &Polynomial, b: &Polynomial| -> Polynomial {
    a.clone() + (-b.clone())
});

impl_op_ex!(- |a: &Polynomial, b: &Term | -> Polynomial {
    a.clone() + (-b.to_pol())
});

impl_op_ex!(- |a: &Term, b: &Polynomial| -> Polynomial {
    a.to_pol() + (-b.clone())
});

// -=
impl_op_ex!(-= |a: &mut Polynomial, b: &Polynomial| {
    for (bk, bv) in &b.terms {
        if let Some(av) = a.terms.get_mut(&bk) {
            *av -= bv;
            if av.is_zero() {
                a.terms.remove(&bk.clone());
            }
        } else {
            a.terms.insert(bk.clone(), - bv.clone());
        }
    }
});

// Neg
impl_op_ex!(- |a: &Polynomial| -> Polynomial {
    let mut pol = Polynomial::new();
    for (k, coef) in &a.terms {
        pol.terms.insert(k.clone(), -coef);
    }
    pol
});

// *
impl_op_ex!(* |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut pol = Polynomial::new();
    for (ik, iv) in &a.terms {
        let i = Term::from(ik, iv);
        for (jk, jv) in &b.terms {
            let j = Term::from(jk, jv);
            let l = i.clone() * j;
            if let Some(lv) = pol.terms.get_mut(&l.monomial) {
                *lv += &l.coef;
            } else {
                pol.terms.insert(l.monomial, l.coef);
            }
        }
    }
    pol
});

impl_op_ex!(* |a: &Polynomial, b: &Term| -> Polynomial {
    a * b.to_pol()
});

impl_op_ex!(* |a: &Term, b: &Polynomial| -> Polynomial {
    a.to_pol() * b
});

// *=
impl_op_ex!(*= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() * b;
    a.terms.clear();
    a.terms = c.terms; 
});

// /
impl_op_ex!(/ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    if b.is_zero() {
        panic!("b.is_zero()");
    }
    else if b.terms.len() == 1 {
        let mut pol = Polynomial::new();
        for (ak, av) in &a.terms {
            let i = Term::from(ak, av);
            let (k, v) = b.terms.iter().next().unwrap();
            let u2 = Term::from(k, v); 
            let u = i / &u2;
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    } else {
        panic!("b.terms.len() >= 2");
    }
});

impl_op_ex!(/ |a: &Polynomial, b: &Term| -> Polynomial {
    a / b.to_pol()
});

impl_op_ex!(/ |a: &Term, b: &Polynomial| -> Polynomial {
    a.to_pol() / b
});

// /=
impl_op_ex!(/= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() / b;
    a.terms.clear();
    a.terms = c.terms;
});

impl_op_ex!(/= |a: &mut Polynomial, b: &Term| {
    let c = a.clone() / b.to_pol();
    a.terms.clear();
    a.terms = c.terms;
});

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let  s = self.clone();
        if s.is_zero() {
            return write!(f, "0");
        }
        let mut st = String::new();
        let mut i = 0;
        for (k, coef) in s.terms.iter().rev() {
            if coef > &BigInt::from(0) && i != 0 {
                st.push_str("+ ");
            }
            st.push_str(&Term::from(k, coef).to_string());
            st.push_str(" ");
            i = i + 1;
        }
        write!(f, "{}", st.trim_end())
    }
}

impl<'a> Power<&'a BigInt> for Polynomial {
    fn power(&self, n: &BigInt) -> Self {
        assert!(n >= &Zero::zero(), "n:{}", n.to_string());
        if n == &Zero::zero() {
            TermBuilder::new().coef(1).build().to_pol();
        }
        
        let mut m = self.clone();
        for _i in num_iter::range(BigInt::from(0), n-1) {
            m *= self;
        }
        m
    }
}

impl Power<i64> for Polynomial {
    fn power(&self, n: i64) -> Self {
        let n = BigInt::from(n);
        self.power(&n)
    }
}

impl Polynomial {
    pub fn new() -> Self {
        Polynomial {
            terms: BTreeMap::new(),  
        }
    }

    pub fn square(&self) -> Self {
        self.power(&BigInt::from(2))
    }

    pub fn cube(&self) -> Self {
        self.power(&BigInt::from(3))
    }

    pub fn power_modulo(&self, n: &BigInt, p: &BigInt) -> Self {
        assert!(*n >= Zero::zero());
        let mut b = self.modulo(p);
        let mut r = TermBuilder::new().coef(1).build().to_pol();
        let mut e = n.clone();
        while &e > &One::one() {
            if &e % 2 != Zero::zero() {
                r *= &b;
                r.modular_assign(p);
            }
            b *= b.clone();
            b.modular_assign(p);
            e /= 2;
        }
        r *= b;
        r.modular_assign(p);
        r
    }

    pub fn polynomial_modular(&self, other: &Polynomial, p: &BigInt) -> Self {
        assert!(!other.has_y(), "!other.has_y()");
        let oh = other.highest_term_x();
        let mut r = self.clone();
        loop {
            let rh = r.highest_term_x();
            if rh.xpow() < oh.xpow() {
                break;
            }
            let mut q = TermBuilder::new()
                    .coef(&(rh.coef.clone() * oh.coef.inverse(p)))
                    .xpow(&(rh.xpow() - oh.xpow()))
                    .ypow(&(rh.ypow().clone()))
                    .build();
            q.modular_assign(&p);
            let mut d = q * other;
            d.modular_assign(&p);
            r -= d;
            r.modular_assign(&p);
        }
        r.modular_assign(&p);
        r
    }

    pub fn modulo(&self, p: &BigInt) -> Self {
        let tmp = self.clone();
        let mut pol = Polynomial::new();  
        for (k, coef) in tmp.terms {
            let c = coef.mod_floor(p);
            if p != &Zero::zero() {
                assert!(&c < p, "c {}", &c);
            }
            if !c.is_zero() {
                pol.terms.insert(k, c);
            }
        }

        for (_, coef) in &pol.terms {
            if p != &Zero::zero() {
                assert!(coef < p, "{} {}", &coef, pol.to_string());
            }
        }
        pol
    }

    pub fn modular_assign(&mut self, p: &BigInt) {
        let mut del: BTreeSet<term::Monomial> = BTreeSet::new();
        for (k, v) in &mut self.terms {
            let c = v.mod_floor(p);
            if p != &Zero::zero() {
                assert!(&c < p, "c {}", &c);
            }
            *v = c.clone();
            if c.is_zero() {
                del.insert(k.clone());
            }
        }

        // remove
        for k in &del {
            if let Some(_v) = self.terms.get_mut(&k) {
                self.terms.remove(k);
            }
        }

        // check
        for (_, coef) in &self.terms {
            if p != &Zero::zero() {
                assert!(coef < p, "{} {}", &coef, self.to_string());
            }
        }
    }

    pub fn highest_term_x(&self) -> Term {
        if self.is_zero() {
            return Term::new();
        }
        for (k, v) in self.terms.iter().rev() {
            return Term::from(k, v);
        }
        panic!("highest_term_x assert!");
    }

    pub fn is_zero(&self) -> bool {
        self.terms.len() == 0
    }

    pub fn is_gcd_one(&self, other: &Self, p: &BigInt) -> bool {
        let s = self.highest_term_x();
        let o = other.highest_term_x();
        if s.xpow() < o.xpow() {
            let m = other.polynomial_modular(self, p);
            return !m.is_zero();
        }
        let m = self.polynomial_modular(other, p);
        return !m.is_zero();
    }

    pub fn has_y(&self) -> bool {
        for (k, v) in &self.terms {
            let i = Term::from(k, v);
            if i.has_y() {
                return true;
            }
        }
        false
    }

    pub fn to_frob(&self, n: &BigInt) -> Self {
        let mut pol = Polynomial::new();
        for (k, v) in &self.terms {
            let u = Term::from(k, v);
            let u = u.to_frob(n);
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    }

    pub fn to_y_power(&self, n: &BigInt) -> Self {
        let mut pol = Polynomial::new();
        for (k, v) in &self.terms {
            let u = Term::from(k, v);
            let u = u.to_y_power(n);
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    }

    pub fn reduction(&self, a: &BigInt, b: &BigInt) -> Self {
        let mut t = Polynomial::new();
        for (k, v) in &self.terms {
            let u = Term::from(k, v);
            if u.ypow() >= &BigInt::from(2) {
                let yy = u.ypow().clone().div_floor(&BigInt::from(2));
                let mut e = u.clone().to_pol();
                e /= TermBuilder::new().coef(1).ypow(&(yy.clone() * 2)).build();

                let ee = TermBuilder::new().coef(1).xpow(3).build()
                       + TermBuilder::new().coef(&a.clone()).xpow(1).build()
                       + TermBuilder::new().coef(&b.clone()).build();
                // power() is faster than power_modulo()
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
        for (k, v) in &self.terms {
            let u = Term::from(k, v);
            if u.ypow() >= &BigInt::from(2) {
                let yy = u.ypow().clone().div_floor(&BigInt::from(2));
                let mut e = u.clone().to_pol();
                e /= TermBuilder::new().coef(1).ypow(&(yy.clone() * 2)).build();

                let ee = TermBuilder::new().coef(1).xpow(3).build()
                       + TermBuilder::new().coef(&a.clone()).xpow(1).build()
                       + TermBuilder::new().coef(&b.clone()).build();
                // power() is faster than power_modulo()
                let ee = ee.power_modulo(&yy.clone(), p);
                let ee = ee.reduction_modular(a, b, p);
                e *= ee;
                t += e;
            } else {
                t += u.clone();
            }
        }
        t.modular_assign(p);
        t
    }
}

#[test]
fn polynmomial_test() { 

    use super::termbuilder;
    type TermBuilder = termbuilder::TermBuilder;

    let u1 = TermBuilder::new().coef(1).xpow(4).ypow(2).build();
    let u2 = TermBuilder::new().coef(3).xpow(2).ypow(4).build();
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

    let u21 = u2.power(3);
    assert_eq_str!(u21, "27 x^6 y^12");

    // Modulo
    let u33 = TermBuilder::new().coef(37).xpow(2).ypow(4).build();

    let u35 = TermBuilder::new().coef(35).xpow(3).ypow(5).build();

    let u34 = u33.modulo(&BigInt::from(24));
    assert_eq_str!(u34, "13 x^2 y^4");

    let p36 = u33 + u35;
    let p37 = p36.modulo(&BigInt::from(6));
    assert_eq_str!(p37, "5 x^3 y^5 + x^2 y^4");

    let p38 = p36.modulo(&BigInt::from(5));
    assert_eq_str!(p38, "2 x^2 y^4");

    //println!("{}", line!());

    // power
    let p39 = p38.power(2);
    assert_eq_str!(p39, "4 x^4 y^8");

    let p3_39 = p38.power(3);
    assert_eq_str!(p3_39, "8 x^6 y^12");

    let p40 = p37.power(2);
    assert_eq_str!(p40, "25 x^6 y^10 + 10 x^5 y^9 + x^4 y^8");

    // reduction
    let u41 = TermBuilder::new().coef(1).xpow(1).ypow(3).build();
    let p41 = u41.to_pol();
    let p42 = p41.reduction(&BigInt::from(2), &BigInt::from(7));
    assert_eq_str!(p42, "x^4 y + 2 x^2 y + 7 x y");

    let p422 = TermBuilder::new().coef(2).xpow(3).build().to_pol();
    assert_eq_str!(p422.highest_term_x(), "2 x^3");

    // rem
    let u43 = TermBuilder::new().coef(3).xpow(4).build();
    let u44 = TermBuilder::new().coef(2).xpow(1).build();
    let p47 = u43 + u44;

    let u45 = TermBuilder::new().coef(2).xpow(2).build();
    let u46 = TermBuilder::new().coef(3).build();
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
    let mut p80 = TermBuilder::new().coef(6).xpow(6).build()
                + TermBuilder::new().coef(12).xpow(4).build();
    let p81 = TermBuilder::new().coef(2).xpow(6).build().to_pol();
    p80 -= p81;
    assert_eq_str!(p80, "4 x^6 + 12 x^4");

    // /=
    let u50 = TermBuilder::new().coef(6).xpow(6).build();
    let u51 = TermBuilder::new().coef(12).xpow(4).build();
    let u52 = TermBuilder::new().coef(3).xpow(2).build();
    let mut p50 = &u50 + &u51;
    let p51 = u52.to_pol();
    p50 /= &p51;
    assert_eq_str!(p50, "2 x^4 + 4 x^2");
}

#[test]
fn isogeny_test() {
    use super::termbuilder;
    type TermBuilder = termbuilder::TermBuilder;

    let a = BigInt::from(1132);
    let b = BigInt::from(278);
    let pp = BigInt::from(2003);

    let e1 = TermBuilder::new().coef(1).ypow(2).build();

    let p = TermBuilder::new().coef(1).xpow(2).build()
    + TermBuilder::new().coef(301).xpow(1).build()
    + TermBuilder::new().coef(527).build();

    let q = TermBuilder::new().coef(1).xpow(1).build()
    + TermBuilder::new().coef(301).build();

    let r = TermBuilder::new().coef(1).xpow(2).build()
    + TermBuilder::new().coef(602).xpow(1).build()
    + TermBuilder::new().coef(1942).build();

    let s = TermBuilder::new().coef(1).xpow(2).build()
    + TermBuilder::new().coef(602).xpow(1).build()
    + TermBuilder::new().coef(466).build();

    let sum = &e1 * r.power(2) * q.power(3) 
    - p.power(3) * s.power(2)
    - TermBuilder::new().coef(500).build() * p * s.power(2) * q.power(2)
    - TermBuilder::new().coef(1005).build() * s.power(2) * q.power(3);
    let sum = sum.reduction_modular(&a, &b, &pp);
    assert_eq_str!(sum, "0");
}
