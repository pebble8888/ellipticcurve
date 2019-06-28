extern crate num_bigint;
extern crate num_traits;

use num_bigint::BigInt;
use num_traits::One;
use std::cmp::Ordering;
use std::ops::Add;

#[derive(Debug, Clone)]
struct Polynomial<'a> {
    units: Vec<&'a Unit<'a>>,
}

// 3x^4y^5
// coef: 3
// variables: [x^4, y^5]
#[derive(Debug, Clone)]
struct Unit<'a> {
    coef: BigInt,
    variables: Vec<&'a Variable<'a>>, // always sorted
}

// x^4
// power: 4
// name: x
#[derive(Debug, Clone, PartialEq, Eq)]
struct Variable<'a> {
    power: BigInt,
    name: &'a str,
}

use std::fmt;

// Unit

/*
impl<'a> Add for Unit<'a> {
    type Output = Unit
    fn add(self, other: Output) -> Output {

    }
}
*/

impl<'a> PartialEq for Unit<'a> {
    fn eq(&self, other: &Self) -> bool {
        let mut s = self.clone();
        s.normalize();
        let mut o = other.clone();
        o.normalize();
        s.coef == o.coef && s.variables == s.variables
    }
}

// TODO:
/*
impl<'a> Ord for Unit<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
    }
}
*/

impl<'a> fmt::Display for Polynomial<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = self.clone();
        s.normalize();
        let mut st = String::new();
        for i in &s.units {
            st.push_str(&i.to_string());
            st.push_str(" ");
        }
        write!(f, "{}", st.trim_end())
    }
}

impl<'a> Polynomial<'a> {
    pub fn normalize(&mut self) {
        self.units.sort();
    }
}

impl<'a> Unit<'a> {
    pub fn normalize(&mut self) {
        self.variables.sort();
    }
}

impl<'a> fmt::Display for Unit<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = self.clone();
        s.normalize();
        if s.coef == One::one() {
            let mut st = String::new();
            for i in &s.variables {
                st.push_str(&i.to_string());
                st.push_str(" ");
            }
            write!(f, "{}", st.trim_end())
        } else {
            let mut st = String::new();
            st.push_str(&s.coef.to_string());
            st.push_str(" ");
            for i in &s.variables {
                st.push_str(&i.to_string());
                st.push_str(" ");
            }
            write!(f, "{}", st.trim_end())
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
    assert_eq!(u1.to_string(), "a^2 b^2");
    assert_eq!(u2.to_string(), "a^2 b^2");

    let u3_1 = Unit {
        coef: BigInt::from(3),
        variables: vec![&a2, &b2],
    };
    let u5_2 = Unit {
        coef: BigInt::from(5),
        variables: vec![&b2, &a2],
    };
    assert_eq!(u3_1.to_string(), "3 a^2 b^2");
    let p1 = Polynomial {
        units: vec![&u3_1, &u5_2],
    };
    assert_eq!(p1.to_string(), "");
}
