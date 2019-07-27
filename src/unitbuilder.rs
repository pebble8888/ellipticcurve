
use num_bigint::BigInt;
use std::default::Default;
use super::unit;

type Unit = unit::Unit;

#[derive(Debug, Clone, Default)]
pub struct UnitBuilder {
    coef: i64,
    xpow: i64,
    ypow: i64,
}

impl UnitBuilder {
    pub fn new() -> UnitBuilder {
        Default::default()
    }

    pub fn coef(&mut self, coef: i64) -> &mut UnitBuilder {
        self.coef = coef;
        self
    }

    pub fn xpow(&mut self, xpow: i64) -> &mut UnitBuilder {
        self.xpow = xpow; 
        self
    }

    pub fn ypow(&mut self, ypow: i64) -> &mut UnitBuilder {
        self.ypow = ypow;
        self
    }

    pub fn finalize(&self) -> Unit {
        Unit { 
            coef: BigInt::from(self.coef),
            xpow: BigInt::from(self.xpow),
            ypow: BigInt::from(self.ypow),
        }
    }
}

