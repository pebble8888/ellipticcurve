
use num_bigint::BigInt;
use std::default::Default;
use super::unit;

type Unit = unit::Unit;
type UnitKey = unit::UnitKey;

#[derive(Debug, Clone, Default)]
pub struct UnitBuilder {
    coef: BigInt,
    xpow: BigInt,
    ypow: BigInt,
}

impl UnitBuilder {
    pub fn new() -> UnitBuilder {
        Default::default()
    }

    pub fn coef(&mut self, coef: &BigInt) -> &mut UnitBuilder {
        self.coef = coef.clone();
        self
    }

    pub fn coef_i(&mut self, coef: i64) -> &mut UnitBuilder {
        self.coef = BigInt::from(coef);
        self
    }

    pub fn xpow(&mut self, xpow: &BigInt) -> &mut UnitBuilder {
        self.xpow = xpow.clone();
        self
    }

    pub fn xpow_i(&mut self, xpow: i64) -> &mut UnitBuilder {
        self.xpow = BigInt::from(xpow); 
        self
    }

    pub fn ypow(&mut self, ypow: &BigInt) -> &mut UnitBuilder {
        self.ypow = ypow.clone();
        self
    }

    pub fn ypow_i(&mut self, ypow: i64) -> &mut UnitBuilder {
        self.ypow = BigInt::from(ypow);
        self
    }

    pub fn finalize(&self) -> Unit {
        Unit {
          coef: self.coef.clone(),
          key: UnitKey {
              xpow: self.xpow.clone(),
              ypow: self.ypow.clone(),
          },
        }
    }
}

