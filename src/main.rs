#[allow(dead_code, unused_variables)]
mod constants;
use constants::*;
use prettytable::Table;

#[macro_use]
extern crate prettytable;

// const A1: [f64; 5] = [2.36e+05, 1.90e+05, 9.2e+04, 4.23e+07, 1.49e+07];
const A1: f64 = 5.7718e7;
const A2: f64 = 4.87e7;
const N1: f64 = 1.08e16;
const N2: f64 = 6.52e15;
const RHS: f64 = A1 * N1 / (A2 * N2);
const EPS: f64 = 10.0 * f64::EPSILON;
// const EPS: f64 = 1.0e-4;
const MAX_ITERATIONS: usize = 10000;
const EVS: &[f64] = &[30.0, 40.0, 50.0, 60.0, 70.0];
const OUTPUT_FAILED: bool = false;

fn main() {
    run();
}

fn run() {
    let mut table = Table::new();
    let matrix: Vec<Vec<(LineSpec, LineSpec)>> = create_matrix(EVS.to_vec());

    matrix.iter().enumerate().for_each(|(i, line_specs)| {
        let ev = EVS[i];
        table.add_row(row![format!("{ev:?} eV"), "", ""]);
        line_specs.iter().for_each(|(left, right)| {
            let result = newton(left.1, right.1);
            if let Ok((te, diff)) = result {
                table.add_row(row![
                    format!("CS{:?} - CS{:?}", left.0, right.0),
                    format!("{:?}", te),
                    format!("{:?}", diff)
                ]);
            } else if let Err((te, diff)) = result {
                if OUTPUT_FAILED {
                    table.add_row(row![
                        format!("CS{:?} - CS{:?}", left.0, right.0),
                        format!("{:?} (Failed)", te),
                        format!("{:?}", diff),
                    ]);
                }
            }
        });
        table.add_row(row!["", "", ""]);
    });

    table.printstd();
}

fn newton(cs1: (&[f64], &[f64]), cs2: (&[f64], &[f64])) -> Result<(f64, f64), (f64, f64)> {
    let mut te: f64 = 1.0;
    let mut diff: f64 = 1.0;
    let mut count: usize = 0;
    while diff.abs() > EPS {
        if count >= MAX_ITERATIONS {
            return Err((
                te,
                diff,
                //     format!(
                //     "Newton's method failed to converge. Please check the initial values and try again."
                // ),
            ));
        }
        (diff, te) = intercept(te, cs1, cs2);
        count += 1;
    }
    Ok((te, diff))
}

type LineSpec = (usize, (&'static [f64], &'static [f64]));

fn create_matrix(evs: Vec<f64>) -> Vec<Vec<(LineSpec, LineSpec)>> {
    let list: Vec<Vec<(&'static [f64], &'static [f64])>> =
        evs.iter().map(|&v| extract_cs_by_ev(v)).collect();

    list.iter()
        .map(|cs_list| {
            // 特定のeVにおけるcs_list
            let mut result: Vec<(LineSpec, LineSpec)> = vec![];
            let len = cs_list.len();

            for i in 0..len {
                for j in 0..len {
                    if i == j {
                        continue;
                    }
                    let left: LineSpec = (i + 1, cs_list[i]);
                    let right: LineSpec = (j + 1, cs_list[j]);

                    result.push((left, right))
                }
            }

            result
        })
        .collect::<Vec<_>>()
}

fn extract_cs_by_ev(ev: f64) -> Vec<(&'static [f64], &'static [f64])> {
    vec![CS1, CS2, CS3, CS4, CS5, CS6]
        .iter()
        .map(|&v| {
            let index = v.0.iter().filter(|&x| x < &ev).count();
            (&v.0[0..index], &v.1[0..index])
        })
        .collect()
}

/// returns next te
fn intercept(te: f64, cs1: (&[f64], &[f64]), cs2: (&[f64], &[f64])) -> (f64, f64) {
    let df = diff(te, cs1, cs2);
    (df, te - (EPS / (diff(te + EPS, cs1, cs2) / df - 1.0)))
}

fn diff(te: f64, cs1: (&[f64], &[f64]), cs2: (&[f64], &[f64])) -> f64 {
    let l = lhs(te, cs1, cs2);
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
    fn test_extract() {
        let result = extract_cs_by_ev(20.0);
        println!("{:?}", result);
    }
}
