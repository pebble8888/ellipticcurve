use num_bigint::BigInt;
use std::default::Default;
use super::term;

type Term = term::Term;
type Monomial = term::Monomial;

#[derive(Debug, Clone, Default)]
pub struct TermBuilder {
    coef_: BigInt,
    xpow_: BigInt,
    ypow_: BigInt,
}

pub trait TermBuildable<T> {
    fn coef(&mut self, coef: T) -> &mut Self;
    fn xpow(&mut self, xpow: T) -> &mut Self;
    fn ypow(&mut self, ypow: T) -> &mut Self;
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
}

impl TermBuilder {
    pub fn new() -> TermBuilder {
        Default::default()
    }

    pub fn build(&self) -> Term {
        Term {
            coef: self.coef_.clone(),
            monomial: Monomial {
                xpow: self.xpow_.clone(),
                ypow: self.ypow_.clone(),
            },
        }
    }
}

