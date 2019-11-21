use std::cmp::Ordering;
use std::fmt;

/// subscripted variable
/// if i = 0 and j = 0, omit it
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct SubscriptedVariable {
    pub i: u64,
    pub j: u64,
    pub empty: bool,
}

impl Default for SubscriptedVariable {
    fn default() -> Self {
        Self { i: 0, j: 0, empty: true }
    }
}

impl Ord for SubscriptedVariable {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.empty && other.empty {
            return Ordering::Equal;
        } else if self.empty && !other.empty {
            return Ordering::Less;
        } else if !self.empty && other.empty {
            return Ordering::Greater;
        }

        if self.i < other.i {
            return Ordering::Less;
        } else if self.i > other.i {
            return Ordering::Greater;
        }

        if self.j < other.j {
            return Ordering::Less;
        } else if self.j > other.j {
            return Ordering::Greater;
        }
        return Ordering::Equal;
    }
}

impl PartialOrd for SubscriptedVariable {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for SubscriptedVariable {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.empty {
            write!(f, "")
        } else {
            write!(f, "c_{}_{}", self.i, self.j)
        }
    }
}

#[test]
fn subscripted_variable_test1() {
    use num_bigint::BigInt;
    use super::term_builder;
    use super::term_builder::TermBuildable;

    let t1 = term_builder::TermBuilder::new()
        .coef(&BigInt::from(1))
        .variable_ij(2u64, 3u64)
        .build()
        .to_pol();
    assert_eq_str!(t1, "c_2_3");

    let t2 = term_builder::TermBuilder::new()
        .coef(&BigInt::from(45))
        .variable_ij(7u64, 8u64)
        .build()
        .to_pol();
    assert_eq_str!(t2, "45 c_7_8");

    let t3 = term_builder::TermBuilder::new()
        .coef(&BigInt::from(671))
        .variable_ij(19u64, 21u64)
        .xpow(5)
        .build()
        .to_pol();
    assert_eq_str!(t3, "671 x^5 c_19_21");

    let t4 = t1 + t2 + t3;
    assert_eq_str!(t4, "671 x^5 c_19_21 + 45 c_7_8 + c_2_3");
}

