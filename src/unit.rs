extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use std::cmp::Ordering;
use std::fmt;
use std::ops; 
use super::polynomial;
use super::bigint::{Power, PowerModular, DivFloor, RemFloor};
use super::unitbuilder;

type Polynomial = polynomial::Polynomial;
type UnitBuilder = unitbuilder::UnitBuilder;

// coef x^xpow y^ypow
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Unit {
    pub coef: BigInt,
    pub key: UnitKey,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct UnitKey {
    pub xpow: BigInt,
    pub ypow: BigInt,
}

impl_op_ex!(+ |a: &Unit, b: &Unit| -> Polynomial {
    let mut pol = a.clone().to_pol();
    if let Some(av) = pol.units.get_mut(&b.key) {
        *av += b.coef.clone(); 
    } else {
        pol.units.insert(b.key.clone(), b.coef.clone());
    }
    pol
});

impl_op_ex!(* |a: &Unit, b: &Unit| -> Unit {
    Unit {
        coef: &a.coef * &b.coef,
        key: UnitKey { 
            xpow: a.xpow() + b.xpow(),
            ypow: a.ypow() + b.ypow(),
        },
    }
});

impl_op_ex!(- |a: &Unit, b: &Unit| -> Polynomial {
    a + (-b)
});

impl_op_ex!(/ |a: &Unit, b: &Unit| -> Unit {
    UnitBuilder::new()
        .coef(&a.coef.div_floor(&b.coef))
        .xpow(&(a.xpow() - b.xpow()))
        .ypow(&(a.ypow() - b.ypow())).
        finalize()
});

impl_op_ex!(- |a: &Unit| -> Unit {
    UnitBuilder::new()
        .coef(&(-a.coef.clone()))
        .xpow(&a.xpow().clone())
        .ypow(&a.ypow().clone())
        .finalize()
});

impl Unit {
    pub fn from(key: &UnitKey, coef: &BigInt) -> Self {
        Unit {
            coef: coef.clone(),
            key: key.clone(),
        }
    }

    pub fn xpow(&self) -> &BigInt {
        &self.key.xpow
    }

    pub fn ypow(&self) -> &BigInt {
        &self.key.ypow
    }

    pub fn new() -> Unit {
        UnitBuilder::new().finalize()
    }

    pub fn to_pol(&self) -> Polynomial {
        let mut pol = Polynomial::new();
        let u = self.clone();
        pol.units.insert(u.key, u.coef);
        pol
    }

    pub fn equal_order(&self, other: &Self) -> bool {
        return &self.xpow() == &other.xpow()
            && &self.ypow() == &other.ypow()
    }

    pub fn power(&self, n: &BigInt) -> Self {
        UnitBuilder::new()
          .coef(&self.coef.clone().power(&n.clone()))
          .xpow(&(self.xpow() * n.clone()))
          .ypow(&(self.ypow() * n.clone()))
          .finalize()
    }

    pub fn power_i(&self, n: i64) -> Self {
        self.power(&BigInt::from(n))
    }

    pub fn power_modular(&self, n: &BigInt, p: &BigInt) -> Self {
        UnitBuilder::new()
            .coef(&self.coef.power_modular(n, p))
            .xpow(&(self.xpow() * n.clone()))
            .ypow(&(self.ypow() * n.clone()))
            .finalize()
    }

    pub fn to_frob(&self, n: &BigInt) -> Self {
        UnitBuilder::new()
          .coef(&self.coef.clone())
          .xpow(&(self.xpow() * n.clone()))
          .ypow(&(self.ypow() * n.clone()))
          .finalize()
    }

    pub fn to_y_power(&self, n: &BigInt) -> Self {
        UnitBuilder::new()
          .coef(&self.coef.clone())
          .xpow(&self.xpow().clone())
          .ypow(&(self.ypow() * n.clone()))
          .finalize()
    }

    pub fn modular(&self, n: &BigInt) -> Self {
        if n == &Zero::zero() {
            panic!("modular zero!");
        } else {
            UnitBuilder::new()
              .coef(&self.coef.rem_floor(n))
              .xpow(&self.xpow().clone())
              .ypow(&self.ypow().clone())
              .finalize()
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coef == Zero::zero()
    }

    pub fn has_y(&self) -> bool {
        self.ypow() != &BigInt::from(0)
    }
}

impl Ord for Unit {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

impl Ord for UnitKey {
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

impl PartialOrd for UnitKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialOrd for Unit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Unit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.coef == BigInt::from(1) {
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
    let u0: Unit = Default::default();
    assert_eq_str!(u0, "0");

    let u1 = Unit {
        coef: BigInt::from(3),
        .. Default::default()
    };
    assert_eq_str!(u1, "3");

    let u2 = Unit {
        coef: BigInt::from(4),
        .. Default::default()
    };
    let u3 = Unit {
        coef: BigInt::from(3),
        .. Default::default()
    };
    assert_eq_str!(u2.clone() * u3.clone(), "12");
    assert_eq_str!(&u2 * &u3, "12");
    assert_eq_str!(UnitBuilder::new()
        .coef_i(20)
        .xpow_i(1)
        .finalize()
        .modular(&BigInt::from(19)), "x");

    // power_modular

    let u11 = UnitBuilder::new().coef_i(1).xpow_i(4).ypow_i(2).finalize();
    let u8 = - &u11;
    assert_eq_str!(u8, "- x^4 y^2"); 
}


