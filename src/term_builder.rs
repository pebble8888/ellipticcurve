use num_bigint::BigInt;
use std::default::Default;
use super::term;
use super::subscripted_variable;

type Term = term::Term;
type SubscriptedVariable = subscripted_variable::SubscriptedVariable;

#[derive(Debug, Clone, Default)]
pub struct TermBuilder {
    coef_: BigInt,
    xpow_: i32,
    ypow_: i32,
    qpow_: i32,
    variable_: SubscriptedVariable, 
}

pub trait TermBuildable<T> {
    /// coeficient
    fn coef(&mut self, coef: T) -> &mut Self;
}

impl<'a> TermBuildable<&'a BigInt> for TermBuilder {
    fn coef(&mut self, coef: &BigInt) -> &mut TermBuilder {
        self.coef_ = coef.clone();
        self
    }
}

impl<'a> TermBuildable<BigInt> for TermBuilder {
    fn coef(&mut self, coef: BigInt) -> &mut TermBuilder {
        self.coef_ = coef.clone();
        self
    }
}

impl TermBuildable<i32> for TermBuilder {
    fn coef(&mut self, coef: i32) -> &mut TermBuilder {
        self.coef_ = BigInt::from(coef);
        self
    }
}

impl TermBuilder {
    pub fn new() -> TermBuilder {
        TermBuilder {
            coef_: BigInt::from(1),
            .. Default::default()
        }
    }

    pub fn xpow(&mut self, xpow: i32) -> &mut TermBuilder {
        self.xpow_ = xpow;
        self
    }

    pub fn ypow(&mut self, ypow: i32) -> &mut TermBuilder {
        self.ypow_ = ypow;
        self
    }

    pub fn qpow(&mut self, qpow: i32) -> &mut TermBuilder {
        self.qpow_ = qpow;
        self
    }

    pub fn variable(&mut self, variable: SubscriptedVariable) -> &mut TermBuilder {
        self.variable_ = variable;
        self
    }

    pub fn variable_ij(&mut self, i: i32, j: i32) -> &mut TermBuilder {
        assert!(i >= 0);
        assert!(j >= 0);
        self.variable_.i = i;
        self.variable_.j = j;
        self.variable_.empty = false;
        self
    }

    pub fn build(&self) -> Term {
        Term {
            coef: self.coef_.clone(),
            monomial: term::Monomial {
                xpow: self.xpow_.clone(),
                ypow: self.ypow_.clone(),
                qpow: self.qpow_.clone(),
                variable: self.variable_,
            },
        }
    }
}

