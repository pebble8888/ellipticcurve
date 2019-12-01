use crate::bigint::Power;
use num_bigint::BigInt;
use num_integer::Integer;
use std::fmt;
use super::polynomial;
use super::term_builder::TermBuildable;
use super::term_builder;
use crate::bigint::Inverse;
use num_traits::Zero;
use num_traits::One;

type TermBuilder = term_builder::TermBuilder;

/// y^2 = x^3 + a x + b
/// GF(p)
pub struct EllipticCurve {
    pub a: BigInt,
    pub b: BigInt,
    pub p: BigInt,
    pol: polynomial::Polynomial,
    /// rational points
    pub points: Vec<ECPoint> // usize
}

/// Jacobian coordinates point
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ECPoint {
    pub x: BigInt,
    pub y: BigInt,
    pub z: BigInt,
}

impl EllipticCurve {
    pub fn new(a: &BigInt, b: &BigInt, p: &BigInt) -> EllipticCurve {
        let pol = TermBuilder::new().xpow(3).build()
        + TermBuilder::new().coef(a).xpow(1).build()
        + TermBuilder::new().coef(b).build();
        let mut ec = EllipticCurve {
            a: a.clone(),
            b: b.clone(),
            p: p.clone(),
            pol: pol.clone(), 
            points: Vec::new(),
        };
        ec.create_points();
        ec
    }

    pub fn new_raw(a: &BigInt, b: &BigInt, p: &BigInt) -> EllipticCurve {
        let pol = TermBuilder::new().xpow(3).build()
        + TermBuilder::new().coef(a).xpow(1).build()
        + TermBuilder::new().coef(b).build();
        let ec = EllipticCurve {
            a: a.clone(),
            b: b.clone(),
            p: p.clone(),
            pol: pol.clone(), 
            points: Vec::new(),
        };
        ec
    }

    pub fn j_invariant(&self) -> BigInt {
        let n = BigInt::from(4) * self.a.power(3);
        let d = n.clone() + BigInt::from(27) * self.b.power(2); 
        let j = BigInt::from(1728) * n.clone() * d.clone().inverse(&self.p);
        return j.mod_floor(&self.p);
    }

    pub fn is_on_curve(&self, ecpoint: &ECPoint) -> bool {
        if ecpoint.is_infinity() {
            return true;
        }
        let c = self.pol.eval_xy(&ecpoint.x, &ecpoint.y);
        let cc = c.mod_floor(&self.p);
        let d = ecpoint.y.power(2);
        let dd = d.mod_floor(&self.p);
        let diff = cc - dd;
        let r = diff.mod_floor(&self.p);
        r == Zero::zero()
    }

    pub fn canonicalize(&self, point: &ECPoint) -> ECPoint {
        ECPoint {
            x: point.x.mod_floor(&self.p),
            y: point.y.mod_floor(&self.p), 
            z: point.z.mod_floor(&self.p),
        }
    }

    /// Elliptic curve point addition
    pub fn plus(&self, point1: &ECPoint, point2: &ECPoint) -> ECPoint {
        assert!(self.is_on_curve(point1));
        assert!(self.is_on_curve(point2));
        if point1.is_infinity() {
            return ECPoint::new(&point2.x, &point2.y);
        } else if point2.is_infinity() {
            return ECPoint::new(&point1.x, &point2.y);
        }
        let p1 = self.canonicalize(point1);
        let p2 = self.canonicalize(point2);
        let x1 = p1.x.clone();
        let x2 = p2.x.clone();
        let y1 = p1.y.clone();
        let y2 = p2.y.clone();
        if p1.x != p2.x {
            let m = (y2.clone() - y1.clone()) * (x2.clone() - x1.clone()).inverse(&self.p);
            let x3 = m.power(BigInt::from(2)) - x1.clone() - x2.clone();
            let y3 = m * (x1.clone() - x3.clone()) - y1.clone(); 
            return ECPoint::new(
                &x3.mod_floor(&self.p),
                &y3.mod_floor(&self.p));
        } else {
            if y1 != y2 || y1 == Zero::zero() {
                ECPoint::infinity()
            } else {
                let m = (BigInt::from(3) * x1.power(BigInt::from(2)) + self.a.clone()) * (BigInt::from(2) * y1.clone()).inverse(&self.p);
                let x3 = m.power(BigInt::from(2)) - BigInt::from(2) * x1.clone();
                let y3 = m * (x1.clone() - x3.clone()) - y1.clone();
                return ECPoint::new(
                    &x3.mod_floor(&self.p),
                    &y3.mod_floor(&self.p));
            }
        }
    }

    /// Point negation: -P
    pub fn negate(&self, point: &ECPoint) -> ECPoint {
        if point.is_infinity() {
            return point.clone();
        }
        ECPoint::new(
            &point.x,
            &((-point.clone().y).mod_floor(&self.p)))
    }

    /// create rational points 
    pub fn create_points(&mut self) {
        for x in num_iter::range(BigInt::from(0), self.p.clone()) {
            let xpol = x.clone().power(3)
                + self.a.clone() * x.clone()
                + self.b.clone();
            let xpol = xpol.mod_floor(&self.p);
            for y in num_iter::range(BigInt::from(0), self.p.clone()) {
                let ypol = y.clone().power(2);
                let ypol = ypol.mod_floor(&self.p);
                if ypol == xpol {
                    let point = ECPoint::new(&x, &y);
                    self.points.push(point.clone());
                    let minus_point = self.negate(&point);
                    if minus_point != point {
                        self.points.push(minus_point);
                    }
                    break;
                }
            }
        }
        self.points.push(ECPoint::infinity());
    }

    /// get all rational points
    pub fn points(&self) -> Vec<ECPoint> {
        self.points.clone()
    }

    /// EC order by points count
    pub fn order(&self) -> usize {
        return self.points.len();
    }

    /// n * P 
    pub fn multiply_scalar(&self, point: &ECPoint, n: &BigInt) -> ECPoint {
        if n == &Zero::zero() {
            return ECPoint::infinity();
        } else if n < &Zero::zero() {
            let minus_np = self.multiply_scalar(point, &(-n));
            return self.plus(&ECPoint::infinity(), &minus_np);
        } else if n == &One::one() {
            return point.clone();
        }
        let point_s = self.multiply_scalar(&point, &(n % BigInt::from(2)));
        let point_r = self.multiply_scalar(&point_s, &BigInt::from(2));
        if n.is_odd() {
            return self.plus(&point_r, point);
        } else {
            return point_r;
        }
    }

    /// order of P
    pub fn point_order(&self, point: &ECPoint) -> BigInt {
        if point.is_infinity() {
            return One::one();
        }
        let mut x = BigInt::from(2);
        let mut r_point = point.clone();
        loop {
            r_point = self.plus(point, &r_point);
            if r_point.is_infinity() {
                return x;
            }
            x += 1;
        }
    }
}

impl ECPoint {
    pub fn new(x: &BigInt, y:&BigInt) -> ECPoint {
        ECPoint {
            x: x.clone(),
            y: y.clone(),
            z: One::one(),
        }
    }
    pub fn infinity() -> ECPoint {
        ECPoint {
            x: Zero::zero(),
            y: One::one(),
            z: Zero::zero(),
        }
    }
    pub fn is_infinity(&self) -> bool {
        self.z == Zero::zero()
    }
}

impl fmt::Display for EllipticCurve {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "F_{}: y^2 = {}, order:{}, j:{}",
            self.p, self.pol, self.order(), self.j_invariant()) 
    }
}

impl fmt::Display for ECPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_infinity() {
            write!(f, "O")
        } else {
            write!(f, "({}, {})", self.x, self.y)
        }
    }
}

// heavy
/*
#[test]
fn elliptic_curve_test1() {
    let ec2 = EllipticCurve::new(&BigInt::from(1132), &BigInt::from(278), &BigInt::from(2003));
    assert_eq_str!(ec2, "F_2003: y^2 = x^3 + 1132 x + 278, order:1956, j:171");
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(1120), &BigInt::from(1391))), true);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(1120), &BigInt::from(1392))), false);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(894), &BigInt::from(1425))), true);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(894), &BigInt::from(1426))), false);
    let p2 = ECPoint::new(&BigInt::from(1120), &BigInt::from(1391));
    let q2 = ECPoint::new(&BigInt::from(894), &BigInt::from(1425));
    assert_eq_str!(ec2.plus(&p2, &q2), "(1683, 1388)");
    
    let ec4 = EllipticCurve::new(&BigInt::from(500), &BigInt::from(1005), &BigInt::from(2003));
    assert_eq_str!(ec4, "F_2003: y^2 = x^3 + 500 x + 1005, order:1956, j:515");
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(565), &BigInt::from(302))), true);
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(565), &BigInt::from(303))), false);
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(1818), &BigInt::from(1002))), true);
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(1818), &BigInt::from(1000))), false);
    let p4 = ECPoint::new(&BigInt::from(565), &BigInt::from(302));
    let q4 = ECPoint::new(&BigInt::from(1818), &BigInt::from(1002));
    assert_eq_str!(ec4.plus(&p4, &q4), "(1339, 821)");
}
*/

#[test]
fn elliptic_curve_test2() {
    let ec = EllipticCurve::new(&BigInt::from(1), &BigInt::from(1), &BigInt::from(5));
    assert_eq_str!(ec, "F_5: y^2 = x^3 + x + 1, order:9, j:2");
    let points = ec.points();
    assert_eq!(points.len(), 9);
    assert_eq_str!(points[0], "(0, 1)");
    assert_eq_str!(points[1], "(0, 4)");
    assert_eq_str!(points[2], "(2, 1)");
    assert_eq_str!(points[3], "(2, 4)");
    assert_eq_str!(points[4], "(3, 1)");
    assert_eq_str!(points[5], "(3, 4)");
    assert_eq_str!(points[6], "(4, 2)");
    assert_eq_str!(points[7], "(4, 3)");
    assert_eq_str!(points[8], "O");
}

#[test]
fn elliptic_curve_test3() {
    let ec = EllipticCurve::new(&BigInt::from(1), &BigInt::from(1), &BigInt::from(5));
    assert_eq_str!(ec, "F_5: y^2 = x^3 + x + 1, order:9, j:2");
    assert_eq_str!(ec.point_order(&ec.points[0]), "9");
    assert_eq_str!(ec.point_order(&ec.points[1]), "9");
    assert_eq_str!(ec.point_order(&ec.points[2]), "3");
    assert_eq_str!(ec.point_order(&ec.points[3]), "3");
    assert_eq_str!(ec.point_order(&ec.points[4]), "9");
    assert_eq_str!(ec.point_order(&ec.points[5]), "9");
    assert_eq_str!(ec.point_order(&ec.points[6]), "9");
    assert_eq_str!(ec.point_order(&ec.points[7]), "9");
    assert_eq_str!(ec.point_order(&ec.points[8]), "1");
}

#[test]
fn elliptic_curve_test4() {
    let ec = EllipticCurve::new(&BigInt::from(1), &BigInt::from(1), &BigInt::from(29));
    println!("{}", ec);
    for x in num_iter::range(0, ec.order()) {
        println!("P{} {} order {}", x + 1, ec.points[x], ec.point_order(&ec.points[x]));       
    }
}

#[test]
fn elliptic_curve_test5() {
    use primes::PrimeSet;
    use primes::is_prime;

    let mut pset = PrimeSet::new();
    for p in pset.iter().skip(2).take(10) { 
        let ec = EllipticCurve::new(&BigInt::from(1), &BigInt::from(1), &BigInt::from(p));
        print!("{}", ec);
        if is_prime(ec.order() as u64) {
            print!(" order is prime");
        }
        println!("");
    }
}

