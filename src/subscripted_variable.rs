use std::cmp::Ordering;
use std::fmt;
use primes::is_prime;

/// subscripted variable
/// if i = 0 and j = 0, omit it
#[derive(Debug, Clone, Copy)]
pub struct SubscriptedVariable {
    pub i: i64,
    pub j: i64,
    pub empty: bool,
}

impl SubscriptedVariable {
    pub fn new() -> Self {
        Self { i: 0, j: 0, empty: true }
    }
}

impl Default for SubscriptedVariable {
    fn default() -> Self {
        Self { i: 0, j: 0, empty: true }
    }
}

impl PartialEq for SubscriptedVariable {
    fn eq(&self, other: &Self) -> bool {
        if self.empty != other.empty {
            return false;
        }
        if self.empty && other.empty {
            return true;
        }
        self.i == other.i && self.j == other.j
    }
}
impl Eq for SubscriptedVariable {}

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

#[derive(Debug)]
pub struct SubscriptedVariableConverter {
    pub p: i64,
}

impl SubscriptedVariableConverter {
    pub fn new(p: i64) -> Self {
        assert!(p >= 2);
        if !is_prime(p as u64) {
            panic!("p must be prime!");
        }

        SubscriptedVariableConverter {
            p: p,
        }
    }

    pub fn count(&self) -> u64 {
        let mut index: u64 = 0;
        for i in num_iter::range(0, self.p+1) {
            for _j in num_iter::range(i+1, self.p+1) {
                index += 1;
            }
        } 
        for _i in num_iter::range(0, self.p+1) {
            index += 1;
        }
        index
    }

    pub fn index_from_variable(&self, variable: SubscriptedVariable) -> Option<u64> {
        if variable.empty {
            return None as Option<u64>;
        }

        let mut index: u64 = 0;
        for i in num_iter::range(0, self.p+1) {
            for j in num_iter::range(i+1, self.p+1) {
                if variable.i == i &&
                   variable.j == j {
                    return Some(index);
                }
                index += 1;
            }
        } 
        for i in num_iter::range(0, self.p+1) {
            if variable.i == i &&
               variable.j == i {
                return Some(index);
            }
            index += 1;
        }
        panic!("invalid variable!");
    }

    pub fn variable_from_index(&self, index: u64) -> SubscriptedVariable {
        let mut local_index: i64 = 0;
        for i in num_iter::range(0, self.p+1) {
            for j in num_iter::range(i+1, self.p+1) {
                if local_index == index as i64 {
                    return SubscriptedVariable {
                        i: i,
                        j: j,
                        empty: false,
                    };
                }
                local_index += 1;
            }
        } 
        for i in num_iter::range(0, self.p+1) {
            if local_index == index as i64 {
                return SubscriptedVariable {
                    i: i,
                    j: i,
                    empty: false,
                };
            }
            local_index += 1;
        }
        panic!("invalid index!");
    }
}

#[test]
fn subscripted_variable_test1() {
    use super::term_builder;
    use super::term_builder::TermBuildable;

    let t1 = term_builder::TermBuilder::new()
        .variable_ij(2, 3)
        .build()
        .to_pol();
    assert_eq_str!(t1, "c_2_3");

    let t2 = term_builder::TermBuilder::new()
        .coef(45)
        .variable_ij(7, 8)
        .build()
        .to_pol();
    assert_eq_str!(t2, "45 c_7_8");

    let t3 = term_builder::TermBuilder::new()
        .coef(671)
        .variable_ij(19, 21)
        .xpow(5)
        .build()
        .to_pol();
    assert_eq_str!(t3, "671 x^5 c_19_21");

    let t4 = t1 + t2 + t3;
    assert_eq_str!(t4, "671 x^5 c_19_21 + 45 c_7_8 + c_2_3");
}

#[test]
fn subscripted_converter_test1() {
    let converter = SubscriptedVariableConverter::new(3);
    assert_eq!(converter.count(), 10);

    assert_eq_str!(converter.variable_from_index(0), "c_0_1");
    assert_eq_str!(converter.variable_from_index(1), "c_0_2");
    assert_eq_str!(converter.variable_from_index(2), "c_0_3");
    assert_eq_str!(converter.variable_from_index(3), "c_1_2");
    assert_eq_str!(converter.variable_from_index(4), "c_1_3");
    assert_eq_str!(converter.variable_from_index(5), "c_2_3");
    assert_eq_str!(converter.variable_from_index(6), "c_0_0");
    assert_eq_str!(converter.variable_from_index(7), "c_1_1");
    assert_eq_str!(converter.variable_from_index(8), "c_2_2");
    assert_eq_str!(converter.variable_from_index(9), "c_3_3");

    let v5 = converter.index_from_variable(SubscriptedVariable {
            i: 2,
            j: 3,
            empty: false,
        });
    match v5 {
        Some(n) => {
            assert_eq!(n, 5);
        }
        None => {
            assert!(false);
        }
    }
}

