/*
extern crate num_bigint;
extern crate num_traits;
extern crate num_iter;

use num_bigint::BigInt;
use num_traits::pow::Pow;
*/

//type Pow = num_traits::pow::Pow;

/*
impl num_traits::pow::Pow for BigInt {
    type Output = Self;
    fn pow(self, val: BigInt) -> Self {
        let mut t = BigInt::from(1);
        for _i in num_iter::range(BigInt::from(0), val-1) {
            t = &t * self;
        }
        t
    }
}
*/
