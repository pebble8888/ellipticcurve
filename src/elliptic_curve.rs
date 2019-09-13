use crate::bigint::Power;
use num_bigint::BigInt;
use num_integer::Integer;
use std::fmt;
use super::polynomial;
use super::termbuilder::TermBuildable;
use super::termbuilder;
use crate::bigint::Inverse;

type TermBuilder = termbuilder::TermBuilder;

pub struct EllipticCurve {
    pub a: BigInt,
    pub b: BigInt,
    pub p: BigInt,
    pol: polynomial::Polynomial,
}

pub struct ECPoint {
    pub x: BigInt,
    pub y: BigInt,
    pub z: BigInt,
}

impl EllipticCurve {
    pub fn new(a: &BigInt, b: &BigInt, p: &BigInt) -> EllipticCurve {
        let pol = TermBuilder::new().coef(1).xpow(3).build()
        + TermBuilder::new().coef(a).xpow(1).build()
        + TermBuilder::new().coef(b).build();
        EllipticCurve {
            a: a.clone(),
            b: b.clone(),
            p: p.clone(),
            pol: pol.clone(), 
        }
    }

    pub fn is_on_curve(&self, ecpoint: &ECPoint) -> bool {
        if ecpoint.is_infinity() {
            return true;
        }
        let c = self.pol.eval(&ecpoint.x, &ecpoint.y);
        let cc = c.mod_floor(&self.p);
        let d = ecpoint.y.power(2);
        let dd = d.mod_floor(&self.p);
        let diff = cc - dd;
        let r = diff.mod_floor(&self.p);
        r == BigInt::from(0)
    }

    pub fn canonicalize(&self, point: &ECPoint) -> ECPoint {
        ECPoint {
            x: point.x.mod_floor(&self.p),
            y: point.y.mod_floor(&self.p), 
            z: point.z.mod_floor(&self.p),
        }
    }

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
            let x3 = m.power(&BigInt::from(2)) - x1.clone() - x2.clone();
            let y3 = m * (x1.clone() - x3.clone()) - y1.clone(); 
            ECPoint {
                x: x3.mod_floor(&self.p),
                y: y3.mod_floor(&self.p),
                z: BigInt::from(0),
            }
        } else {
            if y1 != y2 || y1 == BigInt::from(0) {
                ECPoint::infinity()
            } else {
                let m = (BigInt::from(3) * x1.power(&BigInt::from(2)) + self.a.clone()) * (BigInt::from(2) * y1.clone()).inverse(&self.p);
                let x3 = m.power(&BigInt::from(2)) - BigInt::from(2) * x1.clone();
                let y3 = m * (x1.clone() - x3.clone()) - y1.clone();
                ECPoint {
                    x: x3.mod_floor(&self.p),
                    y: y3.mod_floor(&self.p),
                    z: BigInt::from(0),
                }
            }
        }
    }
}

impl ECPoint {
    pub fn new(x: &BigInt, y:&BigInt) -> ECPoint {
        ECPoint {
            x: x.clone(),
            y: y.clone(),
            z: BigInt::from(0),
        }
    }
    pub fn infinity() -> ECPoint {
        ECPoint {
            x: BigInt::from(0),
            y: BigInt::from(0),
            z: BigInt::from(1),
        }
    }
    pub fn is_infinity(&self) -> bool {
        self.z != BigInt::from(0)
    }
}

impl fmt::Display for EllipticCurve {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "F_{}: y^2 = {}", self.p, self.pol) 
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

#[test]
fn elliptic_curve_test() {
    let ec2 = EllipticCurve::new(&BigInt::from(1132), &BigInt::from(278), &BigInt::from(2003));
    assert_eq_str!(ec2, "F_2003: y^2 = x^3 + 1132 x + 278");
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(1120), &BigInt::from(1391))), true);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(1120), &BigInt::from(1392))), false);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(894), &BigInt::from(1425))), true);
    assert_eq!(ec2.is_on_curve(&ECPoint::new(&BigInt::from(894), &BigInt::from(1426))), false);
    let p2 = ECPoint::new(&BigInt::from(1120), &BigInt::from(1391));
    let q2 = ECPoint::new(&BigInt::from(894), &BigInt::from(1425));
    assert_eq_str!(ec2.plus(&p2, &q2), "(1683, 1388)");
    
    let ec4 = EllipticCurve::new(&BigInt::from(500), &BigInt::from(1005), &BigInt::from(2003));
    assert_eq_str!(ec4, "F_2003: y^2 = x^3 + 500 x + 1005");
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(565), &BigInt::from(302))), true);
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(565), &BigInt::from(303))), false);
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(1818), &BigInt::from(1002))), true);
    assert_eq!(ec4.is_on_curve(&ECPoint::new(&BigInt::from(1818), &BigInt::from(1000))), false);
    let p4 = ECPoint::new(&BigInt::from(565), &BigInt::from(302));
    let q4 = ECPoint::new(&BigInt::from(1818), &BigInt::from(1002));
    assert_eq_str!(ec4.plus(&p4, &q4), "(1339, 821)");

}
