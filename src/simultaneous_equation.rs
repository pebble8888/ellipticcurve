use num_integer::Integer;
use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;

pub fn solve(matrix: &Vec<Vec<BigInt>>) -> Vec<Vec<BigInt>> {
    let row_count: usize = matrix.len();  
    let col_count: usize = matrix[0].len();
    let mut a: Vec<Vec<BigInt>> = vec![vec![BigInt::from(0); col_count]; row_count];
    for row in num_iter::range(0, row_count) {
        for col in num_iter::range(0, col_count) {
            a[row][col] = matrix[row][col].clone();
        }
    }

    let mut row: Vec<usize> = Vec::new();
    for i in num_iter::range(0, row_count) {
        row.push(i);
    }

    // (i) LU
    for i in num_iter::range(0, col_count - 1) {
        // (I)
        for j in num_iter::range(i, row_count) {
            if a[row[j]][i] != BigInt::from(0) {
                if i != j {
                    // row swap
                    let b = row[i];
                    row[i] = row[j];
                    row[j] = b;
                }
                break;
            }
        }
        let c1 = a[row[i]][i].clone();
        assert!(c1 != BigInt::from(0));
        // (II)
        for j in num_iter::range(i+1, row_count) {
            let c2 = a[row[j]][i].clone();
            if c2 != BigInt::from(0) {
                let lcm = Integer::lcm(&c1, &c2); 
                let cc1 = lcm.clone() / c1.clone();
                //assert_eq!(cc1.clone() * c1.clone(), lcm.clone());
                let cc2 = lcm.clone() / c2.clone();
                //assert_eq!(cc2.clone() * c2.clone(), lcm.clone());
                for k in num_iter::range(0, col_count) {
                    a[row[j]][k] *= cc2.clone();
                }
                for k in num_iter::range(0, col_count) {
                    let val = a[row[i]][k].clone() * cc1.clone();
                    a[row[j]][k] -= val;
                }
            }
        }
    }

    // (ii) Substitute
    let col_prim = col_count -1;
    for i in num_iter::range(1, col_count - 1).rev() {
        // diagonal coef
        let diag = a[row[i]][i].clone();
        let prim = a[row[i]][col_prim].clone();
        //assert_eq!(prim.mod_floor(&diag.clone()), BigInt::from(0));
        if diag != One::one() {
            a[row[i]][i] = 1.into();
            a[row[i]][col_prim] = prim / diag;
        }
        let val = a[row[i]][col_prim].clone();
        for j in num_iter::range(0, i) {
            let t_val = a[row[j]][i].clone();
            a[row[j]][col_prim] -= t_val * val.clone();
            a[row[j]][i] = Zero::zero();
        }
    }
    let mut b: Vec<Vec<BigInt>> = vec![vec![BigInt::from(0); col_count]; row_count];
    for i in num_iter::range(0, row_count) {
        b[i] = a[row[i]].clone();
    }
    b
}
