extern crate ellipticcurve;
use ellipticcurve::modular_polynomial;

fn main() {
    println!("{}", modular_polynomial::modular_polynomial(2));
    println!("{}", modular_polynomial::modular_polynomial(3));
    println!("{}", modular_polynomial::modular_polynomial(5));
}

