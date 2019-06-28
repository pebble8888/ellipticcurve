extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::One;
use std::cmp::Ordering;

#[derive(Debug)]
struct Unit<'a> {
    coef: BigInt, // 整数係数 
    variables: Vec<&'a Variable<'a>>, // 変数リスト(全て乗算したとみなす)
}

#[derive(Debug, PartialEq, Eq)]
struct Variable<'a> {
    power: BigInt, // べき数
    name: &'a str, //変数の文字列
}

use std::fmt;

// Unit

/*
impl<'a> Ord for Unit<'a> {
    fn cmp(&self, other: &Self) -> Ordering {

    }
}
*/

impl<'a> Unit<'a> {
    pub fn normalize(&self) {
        //self.sort()
        // TODO:
    }
}

impl<'a> fmt::Display for Unit<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.normalize();
        if self.coef == One::one() {
            let mut s = String::new();
            for i in &self.variables {
                s.push_str(&i.to_string());
            }
            write!(f, "{}", s)
        } else {
            let mut s = String::new();
            s.push_str(&self.coef.to_string());
            for i in &self.variables {
                s.push_str(&i.to_string());
            }
            write!(f, "{}", s)
        }
    }
}

// Variable

impl<'a> PartialOrd for Variable<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> Ord for Variable<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.name < other.name {
            Ordering::Less
        } 
        else if self.name > other.name {
            Ordering::Greater
        }
        else if self.power < other.power {
            Ordering::Less
        }
        else if self.power > other.power {
            Ordering::Greater
        }
        else { 
            Ordering::Equal
        }
    }
}

impl<'a> fmt::Display for Variable<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.power == One::one() { 
            write!(f, "{}", self.name)
        } else {
            write!(f, "{}^{}", self.name, self.power)
        }
    }
}

// test variable

#[test]
fn variable() {
    let a1 = Variable {
        power: One::one(),
        name: "a",
    };
    let a2 = Variable {
        power: BigInt::from(2),
        name: "a",
    };
    let b1 = Variable {
        power: One::one(),
        name: "b",
    };
    let b2 = Variable {
        power: BigInt::from(2),
        name: "b",
    };
    let a31 = Variable {
        power: BigInt::from(31),
        name: "a",
    };
    let ba3 = Variable {
        power: BigInt::from(3),
        name: "ba",
    };
    assert_eq!(a1.to_string(), "a");
    assert_eq!(a2.to_string(), "a^2");
    assert_eq!(b1.to_string(), "b");
    assert_eq!(b2.to_string(), "b^2");
    assert_eq!(a31.to_string(), "a^31");
    assert_eq!(ba3.to_string(), "ba^3");

    assert!(a1 < a2);
    assert!(a1 < b1);
    assert!(b1 < b2);
}

// test unit

#[test]
fn unit() {
    let a2 = Variable {
        power: BigInt::from(2),
        name: "a",
    };
    let b2 = Variable {
        power: BigInt::from(2),
        name: "b",
    };
    let u1 = Unit {
        coef: One::one(),
        variables: vec![&a2, &b2],
    };
    let u2 = Unit {
        coef: One::one(),
        variables: vec![&b2, &a2],
    };
    assert_eq!(u1.to_string(), "a^2b^2");
    assert_eq!(u2.to_string(), "a^2b^2");
}
