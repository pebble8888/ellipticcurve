use num_integer::Integer;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use std::cmp::Ordering;
use std::fmt;
use std::ops; 
use super::polynomial;
use super::bigint::{Power, PowerModulo};
use super::term_builder;
use super::term_builder::TermBuildable;
use super::subscripted_variable;

type SubscriptedVariable = subscripted_variable::SubscriptedVariable;

/// coef * monomial
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Term {
    pub coef: BigInt,
    pub monomial: Monomial,
}

/// x^xpow * y^ypow * c_i_j * q^qpow
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Monomial {
    pub xpow: BigInt,
    pub ypow: BigInt,
    pub qpow: BigInt,
    pub variable: SubscriptedVariable,
}

impl One for Term {
    fn one() -> Self {
        term_builder::TermBuilder::new().build()   
    }
}

impl_op_ex!(+ |a: &Term, b: &Term| -> polynomial::Polynomial {
    let mut pol = a.clone().to_pol();
    if let Some(av) = pol.terms.get_mut(&b.monomial) {
        *av += b.coef.clone(); 
    } else {
        pol.terms.insert(b.monomial.clone(), b.coef.clone());
    }
    pol
});

impl_op_ex!(* |a: &Term, b: &Term| -> Term {
    if !a.variable().empty && !b.variable().empty {
        panic!("can't multiply another variable!");
    }
    Term {
        coef: &a.coef * &b.coef,
        monomial: Monomial { 
            xpow: a.xpow() + b.xpow(),
            ypow: a.ypow() + b.ypow(),
            qpow: a.qpow() + b.qpow(),
            variable: if !a.variable().empty { a.variable() } else { b.variable() },
        },
    }
});

impl_op_ex!(- |a: &Term, b: &Term| -> polynomial::Polynomial {
    a + (-b)
});

impl_op_ex!(/ |a: &Term, b: &Term| -> Term {
    if !b.variable().empty {
        panic!("b is not zero!");
    }
    term_builder::TermBuilder::new()
        .coef(a.coef.div_floor(&b.coef))
        .xpow(a.xpow() - b.xpow())
        .ypow(a.ypow() - b.ypow())
        .qpow(a.qpow() - b.qpow())
        .variable(a.variable())
        .build()
});

impl_op_ex!(- |a: &Term| -> Term {
    term_builder::TermBuilder::new()
        .coef(-a.coef.clone())
        .xpow(a.xpow())
        .ypow(a.ypow())
        .qpow(a.qpow())
        .variable(a.variable())
        .build()
});

/// Term ^ n
impl Power<i64> for Term {
    fn power(&self, n: i64) -> Self {
        if !self.variable().empty {
            panic!("variable can't power");
        }
        self.power(BigInt::from(n))
    }
}

impl Power<BigInt> for Term {
    fn power(&self, n: BigInt) -> Self {
        self.power(&n)
    }
}

impl<'a> Power<&'a BigInt> for Term {
    fn power(&self, n: &BigInt) -> Self {
        if !self.variable().empty {
            panic!("can't power non zero variable!");
        }
        term_builder::TermBuilder::new()
          .coef(self.coef.power(n))
          .xpow(self.xpow() * n)
          .ypow(self.ypow() * n)
          .qpow(self.qpow() * n)
          .build()
    }
}

impl Monomial {
    pub fn eval_xy(&self, x: &BigInt, y: &BigInt) -> BigInt {
        x.power(&self.xpow) * y.power(&self.ypow)
    }

    pub fn eval_x_polynomial(&self, polynomial: &polynomial::Polynomial) -> polynomial::Polynomial {
        let t = term_builder::TermBuilder::new()
            .ypow(&self.ypow)
            .qpow(&self.qpow)
            .variable(self.variable)
            .build()
            .to_pol();
        polynomial.power(&self.xpow) * t
    }

    pub fn eval_y_polynomial(&self, polynomial: &polynomial::Polynomial) -> polynomial::Polynomial {
        let t = term_builder::TermBuilder::new()
            .xpow(&self.xpow)
            .qpow(&self.qpow)
            .variable(self.variable)
            .build()
            .to_pol();
        polynomial.power(&self.ypow) * t
    }

    pub fn eval_q(&self, q: &BigInt) -> BigInt {
        q.power(&self.qpow)
    }
}

impl Term {
    pub fn from(monomial: &Monomial, coef: &BigInt) -> Self {
        Term {
            coef: coef.clone(),
            monomial: monomial.clone(),
        }
    }

    pub fn xpow(&self) -> &BigInt {
        &self.monomial.xpow
    }

    pub fn ypow(&self) -> &BigInt {
        &self.monomial.ypow
    }

    pub fn qpow(&self) -> &BigInt {
        &self.monomial.qpow
    }

    pub fn variable(&self) -> SubscriptedVariable {
        self.monomial.variable
    }

    pub fn new() -> Term {
        Term {
            coef: Zero::zero(),
            monomial: Monomial {
                xpow: Zero::zero(),
                ypow: Zero::zero(),
                qpow: Zero::zero(),
                variable: SubscriptedVariable {
                    i: 0,
                    j: 0,
                    empty: true,
                },
            },
        }
    }

    /// convert Term to Polynomial
    pub fn to_pol(&self) -> polynomial::Polynomial {
        let mut pol = polynomial::Polynomial::new();
        if !self.coef.is_zero() {
            let u = self.clone();
            pol.terms.insert(u.monomial, u.coef);
        }
        pol
    }

    pub fn is_equal_order(&self, other: &Self) -> bool {
        return &self.xpow() == &other.xpow()
            && &self.ypow() == &other.ypow()
            && &self.qpow() == &other.qpow()
    }

    /// Term^n (mod p)
    pub fn power_modulo(&self, n: &BigInt, p: &BigInt) -> Self {
        term_builder::TermBuilder::new()
            .coef(self.coef.power_modulo(n, p))
            .xpow(self.xpow() * n.clone())
            .ypow(self.ypow() * n.clone())
            .qpow(self.qpow() * n.clone())
            .variable(self.variable())
            .build()
    }

    /// Frobenius map of Term
    /// x -> x^n  y -> y^n
    pub fn to_frob(&self, n: &BigInt) -> Self {
        term_builder::TermBuilder::new()
          .coef(self.coef.clone())
          .xpow(self.xpow() * n.clone())
          .ypow(self.ypow() * n.clone())
          .qpow(self.qpow().clone())
          .variable(self.variable())
          .build()
    }

    /// y -> y^n
    pub fn to_y_power(&self, n: &BigInt) -> Self {
        term_builder::TermBuilder::new()
          .coef(self.coef.clone())
          .xpow(self.xpow().clone())
          .ypow(self.ypow() * n.clone())
          .qpow(self.qpow().clone())
          .variable(self.variable())
          .build()
    }

    /// q -> q^n
    pub fn to_q_power(&self, n: i64) -> Self {
        term_builder::TermBuilder::new()
          .coef(self.coef.clone())
          .xpow(self.xpow().clone())
          .ypow(self.ypow().clone())
          .qpow(self.qpow() * n)
          .variable(self.variable())
          .build()
    }

    /// Term (mod n)
    pub fn modulo(&self, n: &BigInt) -> Self {
        if n.is_zero() {
            panic!("modulo zero!");
        } else {
            term_builder::TermBuilder::new()
              .coef(self.coef.mod_floor(n))
              .xpow(self.xpow().clone())
              .ypow(self.ypow().clone())
              .qpow(self.qpow().clone())
              .variable(self.variable())
              .build()
        }
    }

    /// Term = Term (mod n)
    pub fn modular_assign(&mut self, n: &BigInt) {
        if n.is_zero() {
            panic!("modulo zero!");
        } else {
            self.coef = self.coef.mod_floor(n);
        }
    }

    pub fn has_x(&self) -> bool {
        !self.xpow().is_zero()
    }

    pub fn has_y(&self) -> bool {
        !self.ypow().is_zero()
    }

    pub fn has_q(&self) -> bool {
        !self.qpow().is_zero()
    }

    pub fn eval_x_polynomial(&self, polynomial: &polynomial::Polynomial) -> polynomial::Polynomial {
        let s = term_builder::TermBuilder::new()
            .coef(self.coef.clone())
            .build()
            .to_pol();
        let t = self.monomial.eval_x_polynomial(polynomial);
        s * t
    }

    pub fn eval_y_polynomial(&self, polynomial: &polynomial::Polynomial) -> polynomial::Polynomial {
        let s = term_builder::TermBuilder::new()
            .coef(self.coef.clone())
            .build()
            .to_pol();
        let t = self.monomial.eval_y_polynomial(polynomial);
        s * t
    }

    pub fn derivative_x(&self) -> Term {
        assert!(self.xpow() >= &Zero::zero());
        let mut t = self.clone();
        if t.xpow() > &Zero::zero() { 
            let c = t.xpow().clone(); 
            t.monomial.xpow = c.clone() - BigInt::from(1);
            t.coef *= c;
        } else {
            return term_builder::TermBuilder::new().coef(0).build();
        }
        t
    }

    pub fn derivative_y(&self) -> Term {
        assert!(self.ypow() >= &Zero::zero());
        let mut t = self.clone();
        if t.ypow() > &Zero::zero() { 
            let c = t.ypow().clone(); 
            t.monomial.ypow = c.clone() - BigInt::from(1);
            t.coef *= c;
        } else {
            return term_builder::TermBuilder::new().coef(0).build();
        }
        t
    }
}

impl Ord for Term {
    fn cmp(&self, other: &Self) -> Ordering {
        self.monomial.cmp(&other.monomial)
    }
}

impl Ord for Monomial {
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
        if self.qpow < other.qpow {
            return Ordering::Less;
        } else if self.qpow > other.qpow {
            return Ordering::Greater;
        }
        if self.variable < other.variable {
            return Ordering::Less;
        } else if self.variable > other.variable {
            return Ordering::Greater;
        }
        return Ordering::Equal;
    }
}

impl PartialOrd for Monomial {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialOrd for Term {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Term {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.coef.is_one() {
            if self.xpow().is_zero() &&
               self.ypow().is_zero() &&
               self.qpow().is_zero() {
                if self.variable().empty {
                    write!(f, "1")
                } else {
                    write!(f, "{}", self.variable())
                }
            } else {
                let mut st = String::new();
                if !self.xpow().is_zero() {
                    if self.xpow().is_one() {
                        st.push_str("x");
                    } else {
                        st.push_str("x^");
                        st.push_str(&self.xpow().to_string());
                    }
                    st.push_str(" ");
                }
                if !self.ypow().is_zero() {
                    if self.ypow().is_one() {
                        st.push_str("y");
                    } else {
                        st.push_str("y^");
                        st.push_str(&self.ypow().to_string());
                    }
                    st.push_str(" ");
                }
                if !self.qpow().is_zero() {
                    if self.qpow().is_one() {
                        st.push_str("q");
                    } else {
                        st.push_str("q^");
                        st.push_str(&self.qpow().to_string());
                    }
                    st.push_str(" ");
                }
                if !self.variable().empty {
                    st.push_str(&self.variable().to_string());
                }
                write!(f, "{}", st.trim_end())
            }
        } else if self.coef == BigInt::from(-1) {
            if self.xpow().is_zero() &&
               self.ypow().is_zero() &&
               self.qpow().is_zero() {
                if self.variable().empty {
                    write!(f, "- 1")
                } else {
                    write!(f, "- {}", self.variable())
                }
            } else {
                let mut st = String::new();
                st.push_str("- ");
                if !self.xpow().is_zero() {
                    if self.xpow().is_one() {
                        st.push_str("x");
                    } else {
                        st.push_str("x^");
                        st.push_str(&self.xpow().to_string());
                    }
                    st.push_str(" ");
                }
                if !self.ypow().is_zero() {
                    if self.ypow().is_one() {
                        st.push_str("y");
                    } else {
                        st.push_str("y^");
                        st.push_str(&self.ypow().to_string());
                    }
                    st.push_str(" ");
                }
                if !self.qpow().is_zero() {
                    if self.qpow().is_one() {
                        st.push_str("q");
                    } else {
                        st.push_str("q^");
                        st.push_str(&self.qpow().to_string());
                    }
                    st.push_str(" ");
                }
                if !self.variable().empty {
                    st.push_str(&self.variable().to_string());
                }
                write!(f, "{}", st.trim_end())
            }
        } else {
            let mut st = String::new();
            if self.coef >= Zero::zero() {
                st.push_str(&self.coef.to_string());
            } else {
                let abs_coef = &self.coef * BigInt::from(-1);
                st.push_str("- ");
                st.push_str(&abs_coef.to_string());
            }
            st.push_str(" ");
            if !self.xpow().is_zero() {
                if self.xpow().is_one() {
                    st.push_str("x");
                } else {
                    st.push_str("x^");
                    st.push_str(&self.xpow().to_string());
                }
                st.push_str(" ");
            }
            if !self.ypow().is_zero() {
                if self.ypow().is_one() {
                    st.push_str("y");
                } else {
                    st.push_str("y^");
                    st.push_str(&self.ypow().to_string());
                }
                st.push_str(" ");
            }
            if !self.qpow().is_zero() {
                if self.qpow().is_one() {
                    st.push_str("q");
                } else {
                    st.push_str("q^");
                    st.push_str(&self.qpow().to_string());
                }
                st.push_str(" ");
            }
            if !self.variable().empty {
                st.push_str(&self.variable().to_string());
            }
            write!(f, "{}", st.trim_end())
        }
    }
}

#[test]
fn unit_test() {
    let u0: Term = Default::default();
    assert_eq_str!(u0, "0");

    let u1 = Term {
        coef: BigInt::from(3),
        .. Default::default()
    };
    assert_eq_str!(u1, "3");

    let u2 = Term {
        coef: BigInt::from(4),
        .. Default::default()
    };
    let u3 = Term {
        coef: BigInt::from(3),
        .. Default::default()
    };
    assert_eq_str!(u2.clone() * u3.clone(), "12");
    assert_eq_str!(&u2 * &u3, "12");
    assert_eq_str!(term_builder::TermBuilder::new()
        .coef(20)
        .xpow(1)
        .build()
        .modulo(&BigInt::from(19)), "x");

    // power_modulo

    let u11 = term_builder::TermBuilder::new().coef(1).xpow(4).ypow(2).build();
    let u8 = - &u11;
    assert_eq_str!(u8, "- x^4 y^2"); 
}

