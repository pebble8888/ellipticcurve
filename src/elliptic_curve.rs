use crate::bigint::Power;
use num_bigint::BigInt;
use num_integer::Integer;
use std::fmt;
use std::vec;
use std::ops::Deref;
use super::polynomial;
use super::term_builder::TermBuildable;
use super::term_builder;
use crate::bigint::Inverse;
use num_traits::Zero;
use num_traits::One;
use num_traits::ToPrimitive;

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
        assert!(p >= &BigInt::from(2));
        assert!(primes::is_prime(p.to_u64().unwrap()));
        let pol = term_builder::TermBuilder::new().xpow(3).build()
        + term_builder::TermBuilder::new().coef(a).xpow(1).build()
        + term_builder::TermBuilder::new().coef(b).build();
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
        let pol = term_builder::TermBuilder::new().xpow(3).build()
        + term_builder::TermBuilder::new().coef(a).xpow(1).build()
        + term_builder::TermBuilder::new().coef(b).build();
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
        if !self.is_on_curve(point1) {
            println!("{}", point1);
            assert!(false, "point1 is not on curve");
        }
        assert!(self.is_on_curve(point1), "point1 is not on curve");
        assert!(self.is_on_curve(point2), "point2 is not on curve");
        if point1.is_infinity() {
            return ECPoint::new(&point2.x, &point2.y, &point2.z);
        } else if point2.is_infinity() {
            return ECPoint::new(&point1.x, &point2.y, &point2.z);
        }
        let p1 = self.canonicalize(point1);
        let p2 = self.canonicalize(point2);
        let x1 = p1.x.clone();
        let x2 = p2.x.clone();
        let y1 = p1.y.clone();
        let y2 = p2.y.clone();
        if p1.x != p2.x {
            let m = (&y2 - &y1) * (&x2 - &x1).inverse(&self.p);
            let x3 = m.power(2) - &x1 - &x2;
            let y3 = m * (&x1 - &x3) - &y1; 
            return ECPoint::new(
                &x3.mod_floor(&self.p),
                &y3.mod_floor(&self.p),
                &BigInt::from(1));
        } else {
            if y1 != y2 || y1 == Zero::zero() {
                ECPoint::infinity()
            } else {
                let m = (BigInt::from(3) * x1.power(2) + self.a.clone()) * (BigInt::from(2) * y1.clone()).inverse(&self.p);
                let x3 = m.power(2) - BigInt::from(2) * &x1;
                let y3 = m * (&x1 - &x3) - &y1;
                return ECPoint::new(
                    &x3.mod_floor(&self.p),
                    &y3.mod_floor(&self.p),
                    &BigInt::from(1));
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
            &((-point.clone().y).mod_floor(&self.p)),
            &point.z)
    }

    /// create rational points 
    pub fn create_points(&mut self) {
        for x in num_iter::range(BigInt::from(0), self.p.clone()) {
            let xpol = x.power(3)
                + self.a.clone() * x.clone()
                + self.b.clone();
            let xpol = xpol.mod_floor(&self.p);
            for y in num_iter::range(BigInt::from(0), self.p.clone()) {
                let ypol = y.power(2);
                let ypol = ypol.mod_floor(&self.p);
                if ypol == xpol {
                    let point = ECPoint::new(&x, &y, &BigInt::from(1));
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

    /// EC cardinality by points count
    pub fn cardinality(&self) -> usize {
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
        // TODO:optimise
        let mut pt: ECPoint = point.clone();
        for _ in num_iter::range(BigInt::from(0), n.clone() - BigInt::from(1)) {
            pt = self.plus(&pt, point);
        }
        pt
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

    pub fn division_points(&self, order: &BigInt) -> ECPointVec {
        let mut vec: Vec<ECPoint> = Vec::new();
        for point in &self.points {
            if self.multiply_scalar(point, order) == ECPoint::infinity() { 
                vec.push(point.clone());
            }
        }
        ECPointVec(vec)
    }
}

impl ECPoint {
    pub fn new(x: &BigInt, y:&BigInt, z:&BigInt) -> ECPoint {
        ECPoint {
            x: x.clone(),
            y: y.clone(),
            z: z.clone(),
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
        write!(f, "F_{}: y^2 = {}, cardinality:{}, j:{}",
            self.p, self.pol, self.cardinality(), self.j_invariant()) 
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

pub struct ECPointVec(vec::Vec<ECPoint>);

impl Deref for ECPointVec {
    type Target = vec::Vec<ECPoint>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<ECPointVec> for vec::Vec<ECPoint> {
    fn from(ecpointvec: ECPointVec) -> Self {
        ecpointvec.0.clone()
    }
}

impl From<&ECPointVec> for vec::Vec<ECPoint> {
    fn from(ecpointvec: &ECPointVec) -> Self {
        ecpointvec.0.clone()
    }
}

impl fmt::Display for ECPointVec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "")?;
        let v: vec::Vec<ECPoint> = self.into();
        let mut first = true;
        for i in v {
            if first {
                first = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{}", i)?;
        }
        Ok(())
    }
}

// heavy
/*
#[test]
fn elliptic_curve_test1() {
    let ec2 = EllipticCurve::new(&BigInt::from(1132), &BigInt::from(278), &BigInt::from(2003));
    assert_eq_str!(ec2, "F_2003: y^2 = x^3 + 1132 x + 278, cardinality:1956, j:171");
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(1120), &BigInt::from(1391))), true);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(1120), &BigInt::from(1392))), false);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(894), &BigInt::from(1425))), true);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(894), &BigInt::from(1426))), false);
    let p2 = ECPoint::new(&BigInt::from(1120), &BigInt::from(1391));
    let q2 = ECPoint::new(&BigInt::from(894), &BigInt::from(1425));
    assert_eq_str!(ec2.plus(&p2, &q2), "(1683, 1388)");
    
    let ec4 = EllipticCurve::new(&BigInt::from(500), &BigInt::from(1005), &BigInt::from(2003));
    assert_eq_str!(ec4, "F_2003: y^2 = x^3 + 500 x + 1005, cardinality:1956, j:515");
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
    assert_eq_str!(ec, "F_5: y^2 = x^3 + x + 1, cardinality:9, j:2");
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
    assert_eq_str!(ec, "F_5: y^2 = x^3 + x + 1, cardinality:9, j:2");
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
    for x in num_iter::range(0, ec.cardinality()) {
        println!("P{} {} cardinality {}", x + 1, ec.points[x], ec.point_order(&ec.points[x]));       
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
        if is_prime(ec.cardinality() as u64) {
            print!(" cardinality is prime");
        }
        println!("");
    }
}

#[test]
fn isogeny_test1() {
    let ec = EllipticCurve::new(&BigInt::from(1), &BigInt::from(1), &BigInt::from(19));
    println!("{}", ec);
    for x in num_iter::range(0, ec.cardinality()) {
        let point_order = ec.point_order(&ec.points[x]);
        println!("P{} {} order {}", x + 1, ec.points[x], point_order);
    }
    let points2 = ec.division_points(&BigInt::from(2));
    assert_eq_str!(points2, "O");

    let points3 = ec.division_points(&BigInt::from(3));
    assert_eq_str!(points3, "(2, 7), (2, 12), O");

    let points7 = ec.division_points(&BigInt::from(7));
    assert_eq_str!(points7, "(10, 2), (10, 17), (14, 2), (14, 17), (15, 3), (15, 16), O");
}

#[test]
#[ignore]
fn isogeny_test2() {
    let ec = EllipticCurve::new(&BigInt::from(1811), &BigInt::from(315), &BigInt::from(2003));
    println!("{}", ec);

    let points2 = ec.division_points(&BigInt::from(2));
    assert_eq_str!(points2, "(1164, 0), (1222, 0), (1620, 0), O");

    let points3 = ec.division_points(&BigInt::from(3));
    assert_eq_str!(points3, "(102, 120), (102, 1883), O");

    let points4 = ec.division_points(&BigInt::from(4));
    assert_eq_str!(points4, "(1164, 0), (1222, 0), (1620, 0), O");
}

