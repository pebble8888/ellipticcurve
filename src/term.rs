use num_integer::Integer;
use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use std::cmp::Ordering;
use std::fmt;
use std::ops; 
use super::polynomial;
use super::bigint::{Power, PowerModulo};
use super::term_builder;
use super::term_builder::TermBuildable;

type Polynomial = polynomial::Polynomial;
type TermBuilder = term_builder::TermBuilder;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Term {
    pub coef: BigInt,
    pub monomial: Monomial,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Monomial {
    pub xpow: BigInt,
    pub ypow: BigInt,
}

impl_op_ex!(+ |a: &Term, b: &Term| -> Polynomial {
    let mut pol = a.clone().to_pol();
    if let Some(av) = pol.terms.get_mut(&b.monomial) {
        *av += b.coef.clone(); 
    } else {
        pol.terms.insert(b.monomial.clone(), b.coef.clone());
    }
    pol
});

impl_op_ex!(* |a: &Term, b: &Term| -> Term {
    Term {
        coef: &a.coef * &b.coef,
        monomial: Monomial { 
            xpow: a.xpow() + b.xpow(),
            ypow: a.ypow() + b.ypow(),
        },
    }
});

impl_op_ex!(- |a: &Term, b: &Term| -> Polynomial {
    a + (-b)
});

impl_op_ex!(/ |a: &Term, b: &Term| -> Term {
    TermBuilder::new()
        .coef(&a.coef.div_floor(&b.coef))
        .xpow(&(a.xpow() - b.xpow()))
        .ypow(&(a.ypow() - b.ypow())).
        build()
});

impl_op_ex!(- |a: &Term| -> Term {
    TermBuilder::new()
        .coef(&(-a.coef.clone()))
        .xpow(&a.xpow().clone())
        .ypow(&a.ypow().clone())
        .build()
});

impl Power<i64> for Term {
    fn power(&self, n: i64) -> Self {
        self.power(&BigInt::from(n))
    }
}

impl<'a> Power<&'a BigInt> for Term {
    fn power(&self, n: &BigInt) -> Self {
        TermBuilder::new()
          .coef(&self.coef.clone().power(&n.clone()))
          .xpow(&(self.xpow() * n.clone()))
          .ypow(&(self.ypow() * n.clone()))
          .build()
    }
}

impl Monomial {
    pub fn eval(&self, x:&BigInt, y:&BigInt) -> BigInt {
        x.power(&self.xpow) * y.power(&self.ypow)
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

    pub fn new() -> Term {
        Term {
            coef: Zero::zero(), // BigInt::from(0),
            monomial: Monomial {
                xpow: Zero::zero(), // BigInt::from(0),
                ypow: Zero::zero(), //BigInt::from(0),
            },
        }
    }

    pub fn to_pol(&self) -> Polynomial {
        let mut pol = Polynomial::new();
        let u = self.clone();
        pol.terms.insert(u.monomial, u.coef);
        pol
    }

    pub fn equal_order(&self, other: &Self) -> bool {
        return &self.xpow() == &other.xpow()
            && &self.ypow() == &other.ypow()
    }

    pub fn power_modulo(&self, n: &BigInt, p: &BigInt) -> Self {
        TermBuilder::new()
            .coef(&self.coef.power_modulo(n, p))
            .xpow(&(self.xpow() * n.clone()))
            .ypow(&(self.ypow() * n.clone()))
            .build()
    }

    pub fn to_frob(&self, n: &BigInt) -> Self {
        TermBuilder::new()
          .coef(&self.coef.clone())
          .xpow(&(self.xpow() * n.clone()))
          .ypow(&(self.ypow() * n.clone()))
          .build()
    }

    pub fn to_y_power(&self, n: &BigInt) -> Self {
        TermBuilder::new()
          .coef(&self.coef.clone())
          .xpow(&self.xpow().clone())
          .ypow(&(self.ypow() * n.clone()))
          .build()
    }

    pub fn modulo(&self, n: &BigInt) -> Self {
        if n == &Zero::zero() {
            panic!("modulo zero!");
        } else {
            TermBuilder::new()
              .coef(&self.coef.mod_floor(n))
              .xpow(&self.xpow().clone())
              .ypow(&self.ypow().clone())
              .build()
        }
    }

    pub fn modular_assign(&mut self, n: &BigInt) {
        if n == &Zero::zero() {
            panic!("modulo zero!");
        } else {
            self.coef = self.coef.mod_floor(n);
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coef == Zero::zero()
    }

    pub fn has_y(&self) -> bool {
        self.ypow() != &Zero::zero() // &BigInt::from(0)
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
        if self.coef == One::one() {
            if self.xpow().is_zero() && self.ypow().is_zero() {
                write!(f, "1")
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
                }
                write!(f, "{}", st.trim_end())
            }
        } else if self.coef == BigInt::from(-1) {
            if self.xpow().is_zero() && self.ypow().is_zero() {
                write!(f, "- 1")
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
    assert_eq_str!(TermBuilder::new()
        .coef(20)
        .xpow(1)
        .build()
        .modulo(&BigInt::from(19)), "x");

    // power_modulo

    let u11 = TermBuilder::new().coef(1).xpow(4).ypow(2).build();
    let u8 = - &u11;
    assert_eq_str!(u8, "- x^4 y^2"); 
}


