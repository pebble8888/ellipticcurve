use num_bigint::BigInt;
use std::default::Default;
use super::term;
use super::subscripted_variable;

type Term = term::Term;
type Monomial = term::Monomial;
type SubscriptedVariable = subscripted_variable::SubscriptedVariable;

#[derive(Debug, Clone, Default)]
pub struct TermBuilder {
    coef_: BigInt,
    xpow_: BigInt,
    ypow_: BigInt,
    qpow_: BigInt,
    variable_: SubscriptedVariable, 
}

pub trait TermBuildable<T> {
    /// coeficient
    fn coef(&mut self, coef: T) -> &mut Self;
    /// power of x
    fn xpow(&mut self, xpow: T) -> &mut Self;
    /// power of y
    fn ypow(&mut self, ypow: T) -> &mut Self;
    /// power of q
    fn qpow(&mut self, qpow: T) -> &mut Self;
}

impl<'a> TermBuildable<&'a BigInt> for TermBuilder {
    fn coef(&mut self, coef: &BigInt) -> &mut TermBuilder {
        self.coef_ = coef.clone();
        self
    }
    fn xpow(&mut self, xpow: &BigInt) -> &mut TermBuilder {
        self.xpow_ = xpow.clone();
        self
    }
    fn ypow(&mut self, ypow: &BigInt) -> &mut TermBuilder {
        self.ypow_ = ypow.clone();
        self
    }
    fn qpow(&mut self, qpow: &BigInt) -> &mut TermBuilder {
        self.qpow_ = qpow.clone();
        self
    }
}

impl TermBuildable<i64> for TermBuilder {
    fn coef(&mut self, coef: i64) -> &mut TermBuilder {
        self.coef_ = BigInt::from(coef);
        self
    }
    fn xpow(&mut self, xpow: i64) -> &mut TermBuilder {
        self.xpow_ = BigInt::from(xpow); 
        self
    }
    fn ypow(&mut self, ypow: i64) -> &mut TermBuilder {
        self.ypow_ = BigInt::from(ypow);
        self
    }
    fn qpow(&mut self, qpow: i64) -> &mut TermBuilder {
        self.qpow_ = BigInt::from(qpow);
        self
    }
}

impl TermBuilder {
    pub fn new() -> TermBuilder {
        Default::default()
    }

    pub fn variable(&mut self, variable: SubscriptedVariable) -> &mut TermBuilder {
        self.variable_ = variable;
        self
    }

    pub fn variable_ij(&mut self, i: u64, j: u64) -> &mut TermBuilder {
        self.variable_.i = i;
        self.variable_.j = j;
        self.variable_.empty = false;
        self
    }

    pub fn build(&self) -> Term {
        Term {
            coef: self.coef_.clone(),
            monomial: Monomial {
                xpow: self.xpow_.clone(),
                ypow: self.ypow_.clone(),
                qpow: self.qpow_.clone(),
                variable: self.variable_,
            },
        }
    }
}

