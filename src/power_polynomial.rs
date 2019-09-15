use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Zero, One};
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::{fmt, ops};
use super::bigint::{Inverse, Power};
use super::term;
use super::termbuilder::TermBuildable;
use super::termbuilder;
use super::polynomial;

type Term = term::Term;
type TermBuilder = termbuilder::TermBuilder;
type Polynomial = polynomial::Polynomial;

/*
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PowerPolynomial {
    pub pow: BigInt,
    pub polynomials: BTreeMap<Polynomial>,
}
*/

/*
impl PowerPolynomial {
    pub fn new() -> Self {
        PowerPolynomial {
            values: BTreeMap::new(),
        }
    }
    pub fn from(pol: &Polynomial) -> Self {
        let mut powerpol = PowerPolynomial::new();
        powerpol.values.insert(pol.clone(), BigInt::from(1));
        powerpol
    }
}

// +
impl_op_ex!(+ |a: &PowerPolynomial, b: &PowerPolynomial| -> PowerPolynomial {
    let mut pol = a.clone();
    for (bk, bv) in &b.values {
        if let Some(av) = pol.values.get_mut(&bk) {
            *av += bv;
            if av.is_zero() {
                pol.values.remove(&bk.clone());
            }
        } else {
            pol.values.insert(bk.clone(), bv.clone());
        }
    }
    pol
});

impl_op_ex!(+ |a: &PowerPolynomial, b: &Term| -> PowerPolynomial {
    a + PowerPolynomial::from(&b.to_pol())
});

impl_op_ex!(+ |a: &Term, b: &PowerPolynomial| -> PowerPolynomial {
    PowerPolynomial::from(&a.to_pol()) + b
});

// +=
impl_op_ex!(+= |a: &mut PowerPolynomial, b: &PowerPolynomial| {
    for (bk, bv) in &b.values {
        if let Some(av) = a.values.get_mut(&bk) {
            *av += bv;
            if av.is_zero() {
                a.values.remove(&bk.clone());
            }
        } else {
            a.values.insert(bk.clone(), bv.clone());
        }
    }
});

/*
impl_op_ex!(+= |a: &mut PowerPolynomial, b: &Term| {
    if let Some(av) = a.values.get_mut(&b.monomial) {
        *av += &b.coef;
        if av.is_zero() {
            a.values.remove(&b.monomial.clone());
        }
    } else {
        a.values.insert(b.monomial.clone(), b.coef.clone());
    }
});
*/

// -
impl_op_ex!(- |a: &PowerPolynomial, b: &PowerPolynomial| -> PowerPolynomial {
    a.clone() + (-b.clone())
});

impl_op_ex!(- |a: &PowerPolynomial, b: &Term | -> PowerPolynomial {
    a.clone() + (-PowerPolynomial::from(&b.to_pol()))
});

impl_op_ex!(- |a: &Term, b: &PowerPolynomial| -> PowerPolynomial {
    PowerPolynomial::from(&a.to_pol()) + (-b.clone())
});

// -=
impl_op_ex!(-= |a: &mut PowerPolynomial, b: &PowerPolynomial| {
    for (bk, bv) in &b.values {
        if let Some(av) = a.values.get_mut(&bk) {
            *av -= bv;
            if av.is_zero() {
                a.values.remove(&bk.clone());
            }
        } else {
            a.values.insert(bk.clone(), - bv.clone());
        }
    }
});

/*
// Neg
impl_op_ex!(- |a: &PowerPolynomial| -> PowerPolynomial {
    let mut pol = PowerPolynomial::new();
    for (k, coef) in &a.values {
        pol.values.insert(k.clone(), -coef);
    }
    pol
});
*/
*/
