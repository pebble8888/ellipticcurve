use num_bigint::BigInt;
use std::default::Default;
use super::unit;

type Unit = unit::Unit;
type UnitKey = unit::UnitKey;

#[derive(Debug, Clone, Default)]
pub struct UnitBuilder {
    coef_: BigInt,
    xpow_: BigInt,
    ypow_: BigInt,
}

pub trait UnitBuildable<T> {
    fn coef(&mut self, coef: T) -> &mut Self;
    fn xpow(&mut self, xpow: T) -> &mut Self;
    fn ypow(&mut self, ypow: T) -> &mut Self;
}

impl<'a> UnitBuildable<&'a BigInt> for UnitBuilder {
    fn coef(&mut self, coef: &BigInt) -> &mut UnitBuilder {
        self.coef_ = coef.clone();
        self
    }
    fn xpow(&mut self, xpow: &BigInt) -> &mut UnitBuilder {
        self.xpow_ = xpow.clone();
        self
    }
    fn ypow(&mut self, ypow: &BigInt) -> &mut UnitBuilder {
        self.ypow_ = ypow.clone();
        self
    }
}

impl UnitBuildable<i64> for UnitBuilder {
    fn coef(&mut self, coef: i64) -> &mut UnitBuilder {
        self.coef_ = BigInt::from(coef);
        self
    }
    fn xpow(&mut self, xpow: i64) -> &mut UnitBuilder {
        self.xpow_ = BigInt::from(xpow); 
        self
    }
    fn ypow(&mut self, ypow: i64) -> &mut UnitBuilder {
        self.ypow_ = BigInt::from(ypow);
        self
    }
}

impl UnitBuilder {
    pub fn new() -> UnitBuilder {
        Default::default()
    }

    pub fn build(&self) -> Unit {
        Unit {
            coef: self.coef_.clone(),
            key: UnitKey {
                xpow: self.xpow_.clone(),
                ypow: self.ypow_.clone(),
            },
        }
    }
}

