#[allow(unused_macros)]
macro_rules! assert_eq_str {
    ($left:expr, $right:expr) => ({
        match (&$left, &$right) {
            (left_val, right_val) => {
                let left_val_string = left_val.to_string();
                assert_eq!(&left_val_string, right_val);
            }
        }
    });
}

