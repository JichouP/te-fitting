mod constants;
use constants::*;

// const A1: [f64; 5] = [2.36e+05, 1.90e+05, 9.2e+04, 4.23e+07, 1.49e+07];
const A1: f64 = 5.7718e7;
const A2: f64 = 4.87e7;
const N1: f64 = 1.08e16;
const N2: f64 = 6.52e15;
const RHS: f64 = A1 * N1 / (A2 * N2);
const EPS: f64 = 10.0 * f64::EPSILON;
const MAX_ITERATIONS: usize = 10000;

fn main() {
    let mut te: f64 = 1.0;
    let mut x: f64 = 1.0;
    let mut count: usize = 0;
    while x.abs() > EPS {
        if count >= MAX_ITERATIONS {
            panic!(
                "Newton's method failed to converge. Please check the initial values and try again. {x}, {te}"
            )
        }
        (x, te) = intercept(te);
        count += 1;
    }
    println!("{x}, {te}");
}

/// returns next te
fn intercept(te: f64) -> (f64, f64) {
    let fte = f(te);
    (fte, te - (EPS / (f(te + EPS) / fte - 1.0)))
}

fn f(te: f64) -> f64 {
    let l = lhs(te, CS1, CS2);
    l - RHS
}

fn lhs(te: f64, cs1: (&[f64], &[f64]), cs2: (&[f64], &[f64])) -> f64 {
    c(te, cs1) / c(te, cs2)
}

fn c(te: f64, cs: (&[f64], &[f64])) -> f64 {
    let (energies, omegas) = cs;
    integral(
        energies,
        &energies
            .iter()
            .zip(omegas)
            .map(|(&energy, &omega)| energy * omega * (-energy / te).exp())
            .collect::<Vec<_>>(),
    )
}

fn integral(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() {
        panic!("invalid length");
    };

    let len = x.len();

    let mut result = 0.0;

    for i in 0..len - 1 {
        // 台形
        result += (y[i] + y[i + 1]) * (x[i + 1] - x[i]) / 2.0
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_integral() {
        let actual = integral(&[0.0, 1.0, 2.0], &[0.0, 10.0, 2.0]);
        let expected = 11.0;

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_c() {}
}
