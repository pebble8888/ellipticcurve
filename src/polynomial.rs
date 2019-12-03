use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Zero, One};
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::{fmt, ops};
use super::bigint::{Inverse, Power};
use super::term;
use super::term_builder::TermBuildable;
use super::term_builder;
use super::polynomial;
use super::subscripted_variable;

/// Polynomial
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Polynomial {
    pub terms: BTreeMap<term::Monomial, BigInt>
}

// Polynomial + Polynomial
impl_op_ex!(+ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut pol = a.clone();
    for (bk, bv) in &b.terms {
        if let Some(av) = pol.terms.get_mut(&bk) {
            *av += bv;
            if av.is_zero() {
                pol.terms.remove(&bk);
            }
        } else {
            pol.terms.insert(bk.clone(), bv.clone());
        }
    }
    pol
});

// Polynomial + Term
impl_op_ex!(+ |a: &Polynomial, b: &term::Term| -> Polynomial {
    a + b.to_pol()
});

// Term + Polynomial
impl_op_ex!(+ |a: &term::Term, b: &Polynomial| -> Polynomial {
    a.to_pol() + b
});

// Polynomial += Polynomial
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

// Polynomial += Term
impl_op_ex!(+= |a: &mut Polynomial, b: &term::Term| {
    if let Some(av) = a.terms.get_mut(&b.monomial) {
        *av += &b.coef;
        if av.is_zero() {
            a.terms.remove(&b.monomial.clone());
        }
    } else {
        a.terms.insert(b.monomial.clone(), b.coef.clone());
    }
});

// Polynomial - Polynomial
impl_op_ex!(- |a: &Polynomial, b: &Polynomial| -> Polynomial {
    a.clone() + (-b.clone())
});

// Polynomial - Term
impl_op_ex!(- |a: &Polynomial, b: &term::Term | -> Polynomial {
    a.clone() + (-b.to_pol())
});

// Term - Polynomial
impl_op_ex!(- |a: &term::Term, b: &Polynomial| -> Polynomial {
    a.to_pol() + (-b.clone())
});

// Polynomial -= Polynomial
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

// Negate
impl_op_ex!(- |a: &Polynomial| -> Polynomial {
    let mut pol = Polynomial::new();
    for (m, coef) in &a.terms {
        pol.terms.insert(m.clone(), -coef);
    }
    pol
});

// Polynomial * Polynomial
impl_op_ex!(* |a: &Polynomial, b: &Polynomial| -> Polynomial {
    let mut pol = Polynomial::new();
    for (ik, iv) in &a.terms {
        let i = term::Term::from(ik, iv);
        for (jk, jv) in &b.terms {
            let j = term::Term::from(jk, jv);
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

// Polynomial * Term
impl_op_ex!(* |a: &Polynomial, b: &term::Term| -> Polynomial {
    a * b.to_pol()
});

// Term * Polynomial
impl_op_ex!(* |a: &term::Term, b: &Polynomial| -> Polynomial {
    a.to_pol() * b
});

// Polynomial *= Polynomial
impl_op_ex!(*= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() * b;
    a.terms.clear();
    a.terms = c.terms; 
});

// Polynomial / Polynomial
impl_op_ex!(/ |a: &Polynomial, b: &Polynomial| -> Polynomial {
    if b.is_zero() {
        panic!("b.is_zero()");
    }
    else if b.terms.len() == 1 {
        let mut pol = Polynomial::new();
        for (ak, av) in &a.terms {
            let i = term::Term::from(ak, av);
            let (m, coef) = b.terms.iter().next().unwrap();
            let u2 = term::Term::from(m, coef); 
            let u = i / &u2;
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    } else {
        panic!("b.terms.len() >= 2");
    }
});

// Polynomial / Term
impl_op_ex!(/ |a: &Polynomial, b: &term::Term| -> Polynomial {
    a / b.to_pol()
});

// Term / Polynomial
impl_op_ex!(/ |a: &term::Term, b: &Polynomial| -> Polynomial {
    a.to_pol() / b
});

// Polynomial /= Polynomial
impl_op_ex!(/= |a: &mut Polynomial, b: &Polynomial| {
    let c = a.clone() / b;
    a.terms.clear();
    a.terms = c.terms;
});

// Polynomial /= Term
impl_op_ex!(/= |a: &mut Polynomial, b: &term::Term| {
    let c = a.clone() / b.to_pol();
    a.terms.clear();
    a.terms = c.terms;
});

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = self.clone();
        if s.is_zero() {
            return write!(f, "0");
        }
        let mut st = String::new();
        let mut i = 0;
        for (m, coef) in s.terms.iter().rev() {
            if coef > &BigInt::from(0) && i != 0 {
                st.push_str("+ ");
            }
            if !coef.is_zero() {
                st.push_str(&term::Term::from(m, coef).to_string());
                st.push_str(" ");
            }
            i = i + 1;
        }
        write!(f, "{}", st.trim_end())
    }
}

/// Polynomial ^ BigInt
impl Power<BigInt> for Polynomial {
    fn power(&self, n: BigInt) -> Self {
        self.power(&n)
    }
}
impl<'a> Power<&'a BigInt> for Polynomial {
    fn power(&self, n: &BigInt) -> Self {
        assert!(n >= &Zero::zero(), "n:{}", n.to_string());
        if n.is_zero() {
            return One::one();
        }
        let mut e = n.clone();
        let mut b = self.clone();
        let mut r: Polynomial = One::one();
        while e > One::one() {
            if e.is_odd() {
                r *= b.clone();
            }
            let f = b.clone() * b.clone();
            b = f;
            e /= 2;
        } 
        r * b
    }
}

impl Power<i64> for Polynomial {
    fn power(&self, n: i64) -> Self {
        let n = BigInt::from(n);
        self.power(&n)
    }
}

impl Zero for polynomial::Polynomial {
    fn zero() -> Self {
        polynomial::Polynomial::new()
    }

    fn is_zero(&self) -> bool {
        self.terms.len() == 0
    }
}

impl One for Polynomial {
    fn one() -> Self {
        term::Term::one().to_pol()
    }
}

impl Polynomial {
    pub fn new() -> Self {
        Polynomial {
            terms: BTreeMap::new(),  
        }
    }

    pub fn square(&self) -> Self {
        self.power(2)
    }

    pub fn cube(&self) -> Self {
        self.power(3)
    }

    pub fn power_modulo(&self, n: &BigInt, p: &BigInt) -> Self {
        assert!(*n >= Zero::zero());
        let mut b = self.modulo(p);
        let mut r: Polynomial = One::one();
        let mut e = n.clone();
        while &e > &One::one() {
            if e.is_odd() {
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

    pub fn power_omit_high_order_q(&self, n: i64, order: i64) -> Self {
        let mut pol = polynomial::Polynomial::one();
        // TODO:optimise
        for _i in num_iter::range(0, n) {
            pol *= self.clone();
            pol = pol.omit_high_order_q(order as i64);
        }
        pol
    }

    pub fn polynomial_modular(&self, other: &Polynomial, p: &BigInt) -> Self {
        assert!(!other.has_y(), "!other.has_y()");
        assert!(!other.has_q(), "!other.has_q()");
        let oh = other.highest_term_x();
        let mut r = self.clone();
        loop {
            let rh = r.highest_term_x();
            if rh.xpow() < oh.xpow() {
                break;
            }
            let mut q = term_builder::TermBuilder::new()
                    .coef(&rh.coef * oh.coef.inverse(p))
                    .xpow(rh.xpow() - oh.xpow())
                    .ypow(rh.ypow())
                    .qpow(rh.qpow())
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
        for (m, coef) in tmp.terms {
            let c = coef.mod_floor(p);
            if p != &Zero::zero() {
                assert!(&c < p, "c {}", &c);
            }
            if !c.is_zero() {
                pol.terms.insert(m, c);
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
        for (m, coef) in &mut self.terms {
            let c = coef.mod_floor(p);
            if p != &Zero::zero() {
                assert!(&c < p, "c {}", &c);
            }
            *coef = c.clone();
            if c.is_zero() {
                del.insert(m.clone());
            }
        }

        // remove
        for m in &del {
            if let Some(_v) = self.terms.get_mut(&m) {
                self.terms.remove(m);
            }
        }

        // check
        for (_, coef) in &self.terms {
            if p != &Zero::zero() {
                assert!(coef < p, "{} {}", &coef, self.to_string());
            }
        }
    }

    pub fn highest_term_x(&self) -> term::Term {
        if self.is_zero() {
            return term::Term::new();
        }
        for (m, coef) in self.terms.iter().rev() {
            return term::Term::from(m, coef);
        }
        panic!("highest_term_x assert!");
    }

    pub fn derivative_x(&self) -> Polynomial {
        let mut pol = Polynomial::new();  
        for (m, coef) in &self.terms {
            let t = term::Term::from(&m, &coef).derivative_x();
            if !t.coef.is_zero() {
                pol.terms.insert(t.monomial , t.coef);
            }
        }
        pol
    }

    pub fn derivative_y(&self) -> Polynomial {
        let mut pol = Polynomial::new();  
        for (m, coef) in &self.terms {
            let t = term::Term::from(&m, &coef).derivative_y();
            if !t.coef.is_zero() {
                pol.terms.insert(t.monomial , t.coef);
            }
        }
        pol
    }

    pub fn is_gcd_one(&self, other: &Self, p: &BigInt) -> bool {
        let m = self.gcd(other, p);
        return !m.is_zero();
    }

    pub fn gcd(&self, other: &Self, p: &BigInt) -> Polynomial {
        let s = self.highest_term_x();
        let o = other.highest_term_x();
        if s.xpow() < o.xpow() {
            return other.gcd(self, p);
        }
        let r = self.polynomial_modular(other, p);
        if r.is_zero() {
            return other.to_monic(p);
        }
        return other.gcd(&r, p);
    }

    pub fn to_monic(&self, p: &BigInt) -> Polynomial {
        if self.is_zero() {
            return self.clone();
        }
        let s = self.highest_term_x();
        let inv = s.coef.inverse(p);
        let mut pol = self * term_builder::TermBuilder::new().coef(&inv).build();
        pol.modular_assign(p);
        pol
    }

    pub fn has_x(&self) -> bool {
        for (m, coef) in &self.terms {
            let i = term::Term::from(m, coef);
            if i.has_x() {
                return true;
            }
        }
        false
    }

    pub fn has_y(&self) -> bool {
        for (m, coef) in &self.terms {
            let i = term::Term::from(m, coef);
            if i.has_y() {
                return true;
            }
        }
        false
    }

    pub fn has_q(&self) -> bool {
        for (m, coef) in &self.terms {
            let i = term::Term::from(m, coef);
            if i.has_q() {
                return true;
            }
        }
        false
    }

    /// frobenius map of Polynomial
    /// x -> x^n  y -> y^n
    pub fn to_frob(&self, n: &BigInt) -> Self {
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            let u = term::Term::from(m, coef);
            let u = u.to_frob(n);
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    }

    /// y -> y^n
    pub fn to_y_power(&self, n: &BigInt) -> Self {
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            let u = term::Term::from(m, coef);
            let u = u.to_y_power(n);
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    }

    /// q -> q^n
    pub fn to_q_power(&self, n: i64) -> Self {
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            let u = term::Term::from(m, coef);
            let u = u.to_q_power(n);
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    }

    /// reduction using y^2 = x^3 + a x + b
    pub fn reduction(&self, a: &BigInt, b: &BigInt) -> Self {
        let mut t = Polynomial::new();
        for (m, coef) in &self.terms {
            let u = term::Term::from(m, coef);
            if u.ypow() >= &BigInt::from(2) {
                let yy = u.ypow().clone().div_floor(&BigInt::from(2));
                let mut e = u.clone().to_pol();
                e /= term_builder::TermBuilder::new().ypow(&yy * 2).build();

                let ee = term_builder::TermBuilder::new().xpow(3).build()
                       + term_builder::TermBuilder::new().coef(a.clone()).xpow(1).build()
                       + term_builder::TermBuilder::new().coef(b.clone()).build();
                // power() is faster than power_modulo()
                let ee = ee.power(&yy);
                e *= ee;
                t += e;
            } else {
                t += u.clone();
            }
        }
        t
    }

    /// reduction using y^2 = x^3 + a x + b (mod p)
    pub fn reduction_modular(&self, a: &BigInt, b: &BigInt, p: &BigInt) -> Self {
        let mut t = Polynomial::new();
        for (m, coef) in &self.terms {
            let u = term::Term::from(m, coef);
            if u.ypow() >= &BigInt::from(2) {
                let yy = u.ypow().clone().div_floor(&BigInt::from(2));
                let mut e = u.clone().to_pol();
                e /= term_builder::TermBuilder::new().ypow(yy.clone() * 2).build();

                let ee = term_builder::TermBuilder::new().xpow(3).build()
                       + term_builder::TermBuilder::new().coef(a.clone()).xpow(1).build()
                       + term_builder::TermBuilder::new().coef(b.clone()).build();
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

    /// evaluation using concrete x and y
    pub fn eval_xy(&self, x: &BigInt, y: &BigInt) -> BigInt {
        let mut sum = BigInt::from(0);
        for (m, coef) in &self.terms {
            sum += m.eval_xy(x, y) * coef; 
        }
        sum
    }

    /// evaluation using concrete x
    pub fn eval_x(&self, x: &BigInt) -> Polynomial {
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            let term = term::Term {
                        coef: coef * x.power(&m.xpow),
                        monomial: term::Monomial {
                            xpow: Zero::zero(),
                            ypow: m.ypow.clone(),
                            qpow: m.qpow.clone(),
                            variable: m.variable,
                            }
                       };
            pol += term;
        }
        pol 
    }

    /// evaluation using concrete y
    pub fn eval_y(&self, y: &BigInt) -> Polynomial {
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            let term = term::Term {
                        coef: coef * y.power(&m.ypow),
                        monomial: term::Monomial {
                            xpow: m.xpow.clone(),
                            ypow: Zero::zero(),
                            qpow: m.qpow.clone(),
                            variable: m.variable,
                            }
                       };
            pol += term;
        }
        pol 
    }

    pub fn to_scalar(&self) -> BigInt {
        if self.terms.len() == 0 {
            return Zero::zero();
        }
        if self.terms.len() >= 2 {
            panic!();
        }
        for (m, coef) in &self.terms {
            if !m.xpow.is_zero() {
                panic!();
            }
            if !m.ypow.is_zero() {
                panic!();
            }
            if !m.qpow.is_zero() {
                panic!();
            }
            return coef.clone();
        }
        panic!();
    }

    /// get x degree
    /// assert if y equal not zero or q equal not zero
    pub fn degree_x(&self) -> BigInt {
        assert!(!self.has_y());
        assert!(!self.has_q());
        if self.is_zero() {
            return Zero::zero();
        }
        let s = self.highest_term_x();
        return s.xpow().clone();
    }

    /// omit O(order+1) for q
    pub fn omit_high_order_q(&self, order: i64) -> Polynomial {
        assert!(!self.has_x(), "has_x()");
        assert!(!self.has_y(), "has_y()");
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            if &m.qpow <= &BigInt::from(order) {
                pol.terms.insert(m.clone(), coef.clone());
            }
        }
        pol
    }

    pub fn eval_x_polynomial(&self, polynomial: &polynomial::Polynomial) -> polynomial::Polynomial {
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            let t = term_builder::TermBuilder::new()
                .coef(coef.clone())
                .build()
                .to_pol();
            pol += m.eval_x_polynomial(polynomial) * t;
        }
        pol
    }

    pub fn eval_y_polynomial(&self, polynomial: &polynomial::Polynomial) -> polynomial::Polynomial {
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            let t = term_builder::TermBuilder::new()
                .coef(coef.clone())
                .build()
                .to_pol();
            pol += m.eval_y_polynomial(polynomial) * t;
        }
        pol
    }

    /// get q power coeficient
    pub fn to_q_power_coef(&self, power: &BigInt) -> polynomial::Polynomial {
        assert!(!self.has_x(), "has_x()");
        assert!(!self.has_y(), "has_y()");
        let mut pol = Polynomial::new();
        for (m, coef) in &self.terms {
            if &m.qpow == power {
                let mm = term::Monomial {
                    xpow: Zero::zero(),
                    ypow: Zero::zero(),
                    qpow: Zero::zero(),
                    variable: m.variable,
                };
                pol.terms.insert(mm, coef.clone());
            }
        }
        pol
    }

    /// get variale coeficient
    pub fn to_variable_coef(&self, variable: subscripted_variable::SubscriptedVariable) -> BigInt {
        assert!(!self.has_x(), "has_x()");
        assert!(!self.has_y(), "has_y()");
        assert!(!self.has_q(), "has_q()");
        for (m, coef) in &self.terms {
            if m.variable == variable {
                return coef.clone();
            }
        }
        Zero::zero()
    }
}

#[test]
fn polynmomial_test_zero_one() { 
    assert_eq_str!(Polynomial::one(), "1");
}

#[test]
fn polynmomial_test() { 
    use super::term_builder;
    type TermBuilder = term_builder::TermBuilder;

    let u1 = TermBuilder::new().xpow(4).ypow(2).build();
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
    let u41 = TermBuilder::new().xpow(1).ypow(3).build();
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
fn polynomial_minus_power_test() {
    use super::term_builder;
    type TermBuilder = term_builder::TermBuilder;

    assert_eq_str!(TermBuilder::new().xpow(-4).build().to_pol(), "x^-4");
    let u = TermBuilder::new().coef(3).xpow(-4).build().to_pol() +
            TermBuilder::new().coef(4).xpow(-1).build().to_pol();
    assert_eq_str!(u, "4 x^-1 + 3 x^-4");
}

#[test]
fn polynomial_x_polynomial_test1() {
    use super::term_builder;
    use super::subscripted_variable;

    let a = term_builder::TermBuilder::new().coef(2).qpow(-1).build().to_pol();
    let b = term_builder::TermBuilder::new().coef(3).qpow(2).build().to_pol();
    let qpol = a + b;
    assert_eq_str!(qpol, "3 q^2 + 2 q^-1");

    let m = term::Monomial {
        xpow: BigInt::from(2),
        ypow: BigInt::from(3),
        qpow: BigInt::from(4),
        variable: subscripted_variable::SubscriptedVariable {
            i: 1,
            j: 2,
            empty: false,
        }
    };
    let t = m.eval_x_polynomial(&qpol);
    assert_eq_str!(t, "9 y^3 q^8 c_1_2 + 12 y^3 q^5 c_1_2 + 4 y^3 q^2 c_1_2");
}

#[test]
fn polynomial_y_polynomial_test1() {
    use super::term_builder;
    use super::subscripted_variable;

    let a = term_builder::TermBuilder::new().coef(-3).qpow(-1).build().to_pol();
    let b = term_builder::TermBuilder::new().coef(2).qpow(2).build().to_pol();
    let qpol = a + b;
    assert_eq_str!(qpol, "2 q^2 - 3 q^-1");

    let m = term::Monomial {
        xpow: BigInt::from(2),
        ypow: BigInt::from(3),
        qpow: BigInt::from(4),
        variable: subscripted_variable::SubscriptedVariable {
            i: 1,
            j: 2,
            empty: false,
        }
    };
    let t = m.eval_y_polynomial(&qpol);
    assert_eq_str!(t, "8 x^2 q^10 c_1_2 - 36 x^2 q^7 c_1_2 + 54 x^2 q^4 c_1_2 - 27 x^2 q c_1_2");
}

#[test]
fn polynomial_y_polynomial_test2() {
    use super::term_builder;
    use super::subscripted_variable;

    let a = term_builder::TermBuilder::new().coef(-3).qpow(-1).build().to_pol();
    let b = term_builder::TermBuilder::new().coef(2).qpow(2).build().to_pol();
    let qpol = a + b;
    assert_eq_str!(qpol, "2 q^2 - 3 q^-1");

    let pol = term_builder::TermBuilder::new()
        .coef(2)
        .xpow(2)
        .ypow(3)
        .qpow(4)
        .variable(subscripted_variable::SubscriptedVariable { i: 1, j: 2, empty: false })
        .build()
        + term_builder::TermBuilder::new()
        .coef(3)
        .xpow(1)
        .ypow(1)
        .qpow(4)
        .variable(subscripted_variable::SubscriptedVariable { i: 0, j: 0, empty: false })
        .build();
    let t = pol.eval_y_polynomial(&qpol);
    assert_eq_str!(t, "16 x^2 q^10 c_1_2 - 72 x^2 q^7 c_1_2 + 108 x^2 q^4 c_1_2 - 54 x^2 q c_1_2 + 6 x q^6 c_0_0 - 9 x q^3 c_0_0");
}


#[test]
fn polynomial_zero_power() {
    let qpol = term_builder::TermBuilder::new()
        .qpow(-1)
        .build()
        .to_pol();
    assert_eq_str!(qpol.power(0), "1");
}

#[test]
fn polynomial_x_polynomial_test3() {
    let pol = term_builder::TermBuilder::new()
        .xpow(3)
        .build()
        + term_builder::TermBuilder::new()
        .ypow(3)
        .build();
    assert_eq_str!(pol, "x^3 + y^3");
    let qpol = term_builder::TermBuilder::new()
        .qpow(-1)
        .build()
        .to_pol();
    assert_eq_str!(pol.eval_x_polynomial(&qpol), "y^3 + q^-3");
}

#[test]
fn polynomial_x_polynomial_test4() {
    let pol = term_builder::TermBuilder::new()
        .xpow(3)
        .build()
        + term_builder::TermBuilder::new()
        .ypow(3)
        .build()
        + term_builder::TermBuilder::new()
        .coef(-1)
        .xpow(2)
        .ypow(2)
        .build()
        + term_builder::TermBuilder::new()
        .coef(-1)
        .xpow(1)
        .ypow(1)
        .build();
    assert_eq_str!(pol, "x^3 - x^2 y^2 - x y + y^3");
    let j = term_builder::TermBuilder::new()
        .qpow(-1)
        .build()
        .to_pol();
    let j2 = term_builder::TermBuilder::new()
        .qpow(-2)
        .build()
        .to_pol();
    let p1 = pol.eval_x_polynomial(&j);
    assert_eq_str!(p1, "y^3 - y^2 q^-2 - y q^-1 + q^-3");
    let p2 = p1.eval_y_polynomial(&j2);
    assert_eq_str!(p2, "0");
}

#[test]
fn power_omit_high_order_q_test1() {
    let pol = term_builder::TermBuilder::new().build()
        + term_builder::TermBuilder::new().qpow(1).build();
    assert_eq_str!(pol.power_omit_high_order_q(24, 1), "24 q + 1");
}

#[test]
fn isogeny_test() {
    use super::term_builder;
    type TermBuilder = term_builder::TermBuilder;

    let a = BigInt::from(1132);
    let b = BigInt::from(278);
    let pp = BigInt::from(2003);

    let e1 = TermBuilder::new().ypow(2).build();

    let p = TermBuilder::new().xpow(2).build()
    + TermBuilder::new().coef(301).xpow(1).build()
    + TermBuilder::new().coef(527).build();

    let q = TermBuilder::new().xpow(1).build()
    + TermBuilder::new().coef(301).build();

    let r = TermBuilder::new().xpow(2).build()
    + TermBuilder::new().coef(602).xpow(1).build()
    + TermBuilder::new().coef(1942).build();

    let s = TermBuilder::new().xpow(2).build()
    + TermBuilder::new().coef(602).xpow(1).build()
    + TermBuilder::new().coef(466).build();

    let sum = &e1 * r.power(2) * q.power(3) 
    - p.power(3) * s.power(2)
    - TermBuilder::new().coef(500).build() * p * s.power(2) * q.power(2)
    - TermBuilder::new().coef(1005).build() * s.power(2) * q.power(3);
    let sum = sum.reduction_modular(&a, &b, &pp);
    assert_eq_str!(sum, "0");
}

#[test]
fn derivative_test() {
    use super::term_builder;
    type TermBuilder = term_builder::TermBuilder;
    let p = TermBuilder::new().coef(4).xpow(6).ypow(3).qpow(3).build()
          + TermBuilder::new().coef(2).xpow(5).ypow(2).build();
    assert_eq_str!(p, "4 x^6 y^3 q^3 + 2 x^5 y^2");
    assert_eq_str!(p.derivative_x(), "24 x^5 y^3 q^3 + 10 x^4 y^2");
    assert_eq_str!(p.derivative_y(), "12 x^6 y^2 q^3 + 4 x^5 y");
}

