//! Binary to simulates Quantum Turing Machines (Qtms) and converts them to
//! equivalent quantum circuits

mod circuit;
pub mod cli;
mod complex;
mod gate_g;
mod qtm;
mod register;

use std::{io::stdout, process::ExitCode};

use clap::{CommandFactory, Parser};
use clap_complete::generate;

use crate::cli::{Args, Command};
use crate::complex::C64;
use crate::qtm::{QState, Qtm};
use crate::register::RegQState;

/// Shifts the head right every step without changing state or tape
/// δ(q, a) = |q, a, +1⟩  for all q ∈ Q, a ∈ Γ
fn make_shift_right(num_states: usize, num_symbols: usize) -> Qtm {
  let mut qtm = Qtm::new(num_states, num_symbols);
  for q in 1..=num_states {
    for a in 0..num_symbols {
      qtm.add_transition(q, a, &[(q, a, 1, C64::ONE)]);
    }
  }
  qtm
}

/// Quantum walk on the integers with a two-state Hadamard coin
/// Q = {1,2},  δ(1,a) = (1/√2)(|1,a,+1⟩+|2,a,−1⟩),
///             δ(2,a) = (1/√2)(|1,a,+1⟩−|2,a,−1⟩)
fn make_quantum_walk(num_symbols: usize) -> Qtm {
  let mut qtm = Qtm::new(2, num_symbols);
  let h = C64::INV_SQRT2;
  let mh = -C64::INV_SQRT2;
  for a in 0..num_symbols {
    qtm.add_transition(1, a, &[(1, a, 1, h), (2, a, -1, h)]);
    qtm.add_transition(2, a, &[(1, a, 1, h), (2, a, -1, mh)]);
  }
  qtm
}

/// Classical bit-flip: flips the symbol under the head and moves right
/// Q = {1}, δ(1,0) = |1,1,+1⟩, δ(1,1) = |1,0,+1⟩
fn make_bit_flip() -> Qtm {
  let mut qtm = Qtm::new(1, 2);
  qtm.add_transition(1, 0, &[(1, 1, 1, C64::ONE)]);
  qtm.add_transition(1, 1, &[(1, 0, 1, C64::ONE)]);
  qtm
}

/// Two-state machine: state 1 flips bit and moves right, state 2 passes and
/// moves left
/// δ(1,0)=|2,1,+1⟩, δ(1,1)=|2,0,+1⟩, δ(2,0)=|1,0,−1⟩, δ(2,1)=|1,1,−1⟩
fn make_flip_bounce() -> Qtm {
  let mut qtm = Qtm::new(2, 2);
  qtm.add_transition(1, 0, &[(2, 1, 1, C64::ONE)]);
  qtm.add_transition(1, 1, &[(2, 0, 1, C64::ONE)]);
  qtm.add_transition(2, 0, &[(1, 0, -1, C64::ONE)]);
  qtm.add_transition(2, 1, &[(1, 1, -1, C64::ONE)]);
  qtm
}

/// Prints a full-width section banner with a double-rule border
///
/// Used to visually separate major examples in the output
///
/// # Example output
/// ```text
/// ════════════════════════ … ════════════════════════
///   Example 2 — Quantum walk
/// ════════════════════════ … ════════════════════════
/// ```
fn banner(title: &str) {
  let line = "═".repeat(72);
  println!("\n{line}\n  {title}\n{line}");
}

/// Prints an indented subsection label between two em-dashes
///
/// Used to separate the direct simulation, circuit simulation, and
/// verification blocks within a single example
///
/// # Example output
/// ```text
///   ── Direct simulation ──
/// ```
fn subsection(label: &str) {
  println!("\n  ── {label} ──");
}

/// Prints every step of a direct Qtm simulation history
///
/// For each step, prints the total probability ‖ψ‖² (should stay 1.0 for a
/// valid unitary evolution) followed by every configuration in the
/// superposition with its complex amplitude, formatted via [`QState`]'s
/// `Display` implementation
///
/// # Arguments
/// * `history` — ordered slice of quantum states returned by
///   [`qtm::direct_run`], where `history[0]` is the initial state and
///   `history[t]` is the state after `t` applications of U_δ
fn print_direct_history(history: &[QState]) {
  for (step, state) in history.iter().enumerate() {
    println!("    step {step:>2}  (‖ψ‖² = {:.6})", state.total_prob());
    print!("{state}");
  }
}

/// Prints every step of a circuit simulation history with superposition
/// annotation
///
/// For each step, prints the total probability ‖ψ‖² then, for every distinct
/// active-head location found in the superposition, a human-readable line
/// showing the **integer** tape coordinate, the Qtm state, and the marginal
/// probability of finding the head there
///
/// Register indices live in Z_N = {0, …, N-1}.  This function converts them
/// back to signed integer coordinates centred on 0:
/// * indices 0 …  t   →  +0 … +t  (right half of the tape)
/// * indices t+1 … N-1 →  -(N-t-1) … -1  (left half, wrapped around)
///
/// The raw register state (one [`register::Cell`] per tape square) is then
/// printed via [`RegQState`]'s `Display` implementation
///
/// # Arguments
/// * `history` — ordered slice of register quantum states from
///   [`circuit::circuit_run`]
/// * `t` — number of simulation steps; defines the wrap boundary at index `t`
fn print_circuit_history(history: &[RegQState], t: usize) {
  let n = history
    .last()
    .and_then(|s| s.0.keys().next())
    .map(|rs| rs.n())
    .unwrap_or(0);
  for (step, state) in history.iter().enumerate() {
    println!("    step {step:>2}  (‖ψ‖² = {:.6})", state.total_prob());
    for (idx, qstate, prob) in state.head_distribution() {
      let int_pos: i64 = if idx <= t {
        idx as i64
      } else {
        idx as i64 - n as i64
      };
      println!(
        "           head at tape[{int_pos:+}], Q={qstate}, prob={prob:.6}"
      );
    }
    print!("{state}");
  }
}

/// Verifies that a circuit simulation matches the corresponding direct
/// simulation
///
/// Iterates over every step from 0 to `steps` (inclusive) and delegates to
/// [`circuit::verify_vs_direct`], which compares complex amplitudes
/// configuration-by-configuration.  Tracks the maximum absolute amplitude
/// error across all steps and all basis states, then prints a one-line
/// PASS / FAIL summary
///
/// A result below 1e-8 is considered passing; in practice the error for a
/// correct implementation should be at or below floating-point machine epsilon
/// (~5.5 × 10⁻¹⁶)
///
/// # Arguments
/// * `label`   — short identifier printed in the summary line
/// * `direct`  — step history from [`qtm::direct_run`]
/// * `circuit` — step history from [`circuit::circuit_run`]
/// * `steps`   — total number of steps simulated (length of both slices
///   minus 1)
fn verify_step(
  label: &str,
  direct: &[QState],
  circuit: &[RegQState],
  steps: usize,
) {
  let mut max_err = 0.0_f64;
  for step in 0..=steps {
    let err = circuit::verify_vs_direct(&direct[step], &circuit[step], step);
    if err > max_err {
      max_err = err;
    }
  }
  let ok = if max_err < 1e-8 { "PASS" } else { "FAIL" };
  println!("  {label:15}: max amplitude error = {max_err:.2e}  →  {ok}");
}

fn main() -> ExitCode {
  let args = Args::parse();

  match args.command {
    Some(Command::Completion(a)) => {
      let mut cmd = cli::Args::command();
      let name = cmd.get_name().to_string();
      generate(a.shell, &mut cmd, name, &mut stdout());
      return ExitCode::SUCCESS;
    }
    None => (),
  }

  // ── Example 1: Shift-right ───────────────────────────────────────────────
  banner("Example 1 — Shift-right Qtm  (deterministic, 1 state)");
  {
    let qtm = make_shift_right(1, 2);
    let input = [1usize, 0, 1];
    let steps = 4;
    match qtm.is_valid() {
      Ok(()) => println!("  Bernstein-Vazirani check: PASS"),
      Err(msg) => println!("  INVALID: {msg}"),
    }
    println!("  Tape loop N = {}", circuit::tape_loop_size(steps));

    let direct = qtm::direct_run(&qtm, &input, steps);
    let circ = circuit::circuit_run(&qtm, &input, steps);

    subsection("Direct simulation (U_delta on Config superposition)");
    print_direct_history(&direct);

    subsection("Circuit simulation (Gate-G + F on register state)");
    print_circuit_history(&circ, steps);

    subsection("Amplitude agreement check");
    verify_step("shift-right", &direct, &circ, steps);
  }

  // ── Example 2: Quantum walk ──────────────────────────────────────────────
  banner("Example 2 — Quantum walk with Hadamard coin  (genuinely quantum)");
  {
    let qtm = make_quantum_walk(2);
    let input: [usize; 0] = [];
    let steps = 5;
    match qtm.is_valid() {
      Ok(()) => println!("  Bernstein-Vazirani check: PASS"),
      Err(msg) => println!("  INVALID: {msg}"),
    }
    println!("  Tape loop N = {}", circuit::tape_loop_size(steps));

    let direct = qtm::direct_run(&qtm, &input, steps);
    let circ = circuit::circuit_run(&qtm, &input, steps);

    subsection("Direct simulation");
    print_direct_history(&direct);

    subsection("Circuit simulation");
    print_circuit_history(&circ, steps);

    subsection("Amplitude agreement check");
    verify_step("quantum-walk", &direct, &circ, steps);
  }

  // ── Example 3: Classical bit-flip ────────────────────────────────────────
  banner("Example 3 — Classical bit-flip  (deterministic, 1 state)");
  {
    let qtm = make_bit_flip();
    let input = [1usize, 0, 1];
    let steps = 4;
    match qtm.is_valid() {
      Ok(()) => println!("  Bernstein-Vazirani check: PASS"),
      Err(msg) => println!("  INVALID: {msg}"),
    }
    println!("  Tape loop N = {}", circuit::tape_loop_size(steps));

    let direct = qtm::direct_run(&qtm, &input, steps);
    let circ = circuit::circuit_run(&qtm, &input, steps);

    subsection("Direct simulation");
    print_direct_history(&direct);

    subsection("Circuit simulation");
    print_circuit_history(&circ, steps);

    subsection("Amplitude agreement check");
    verify_step("bit-flip", &direct, &circ, steps);
  }

  // ── Example 4: Flip-bounce ───────────────────────────────────────────────
  banner("Example 4 — Flip-bounce  (2 states: flip+right / pass+left)");
  {
    let qtm = make_flip_bounce();
    let input = [0usize, 1, 0];
    let steps = 4;
    match qtm.is_valid() {
      Ok(()) => println!("  Bernstein-Vazirani check: PASS"),
      Err(msg) => println!("  INVALID: {msg}"),
    }
    println!("  Tape loop N = {}", circuit::tape_loop_size(steps));

    let direct = qtm::direct_run(&qtm, &input, steps);
    let circ = circuit::circuit_run(&qtm, &input, steps);

    subsection("Direct simulation");
    print_direct_history(&direct);

    subsection("Circuit simulation");
    print_circuit_history(&circ, steps);

    subsection("Amplitude agreement check");
    verify_step("flip-bounce", &direct, &circ, steps);
  }

  // ── Complexity table ─────────────────────────────────────────────────────
  banner("Circuit complexity summary  (§6a)");
  println!();
  println!(
    "  t steps on input n ≤ t  →  N = tape_loop_size(t)  (multiple of 3)"
  );
  println!(
    "  Each step: N G-gate applications (depth 3) + N F-flips (depth 1)."
  );
  println!();
  println!(
    "  {:>8}  {:>8}  {:>14}  {:>14}",
    "t", "N", "G-gates/step", "Total G-gates"
  );
  println!("  {}", "-".repeat(52));
  for t in [1, 2, 4, 8, 16, 32, 64, 128] {
    let n = circuit::tape_loop_size(t);
    println!("  {:>8}  {:>8}  {:>14}  {:>14}", t, n, n, n * t);
  }
  println!();
  println!(
    "  Depth  = O(t)   (3 sublayers/step × t steps)      -- improvement over Yao"
  );
  println!(
    "  Size   = O(t²)  (N * t gates ≈ 2t * t)            -- same as Yao"
  );

  ExitCode::SUCCESS
}
