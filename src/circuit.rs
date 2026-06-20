//! Circuit-based simulation (§5, Molina & Watrous 2019)
//!
//! One step = apply G at every centre position (these commute on im(Π))
//!           + apply F₀⋯F_{N-1}

use crate::{
  gate_g::{Local3, apply_g},
  qtm::Qtm,
  register::{RegQState, RegState},
};

/// Smallest multiple of 3 that is ≥ 2t + 1.
pub fn tape_loop_size(t: usize) -> usize {
  let min_n = 2 * t + 1;
  min_n.div_ceil(3) * 3
}

/// Applies gate G centred at one position of the tape loop to the full
/// quantum state
///
/// G is a unitary acting on the three cell-pairs
/// `(S_{centre-1}, T_{centre-1})`, `(S_centre, T_centre)`,
/// `(S_{centre+1}, T_{centre+1})`, with indices taken modulo N. All other
/// cells in the register are left unchanged
///
/// # Arguments
/// * `qtm`    — the Qtm whose transition function defines G.
/// * `qstate` — the current quantum state over register configurations.
/// * `centre` — the tape-loop index (in `0..N`) at which G is centred.
///
/// # Returns
/// The quantum state after G has been applied at `centre`
/// Returns an empty state immediately if `qstate` is empty
fn apply_g_at(qtm: &Qtm, qstate: &RegQState, centre: usize) -> RegQState {
  let n = match qstate.0.keys().next() {
    Some(rs) => rs.n(),
    None => return RegQState::default(),
  };
  let left = (centre + n - 1) % n;
  let right = (centre + 1) % n;

  let mut next = RegQState::default();
  for (rs, &amp) in &qstate.0 {
    let local = Local3([rs.0[left], rs.0[centre], rs.0[right]]);
    for (new_local, g_amp) in apply_g(qtm, local) {
      let mut new_rs = rs.clone();
      new_rs.0[left] = new_local.0[0];
      new_rs.0[centre] = new_local.0[1];
      new_rs.0[right] = new_local.0[2];
      next.add_amp(new_rs, amp * g_amp);
    }
  }
  next.clean();
  next
}

/// Simulates one step of a Qtm computation via the quantum circuit.
///
/// Implements the operator W = (F₀⋯F_{N-1}) V* (F₀⋯F_{N-1}) V, which
/// agrees with AUA* on im(A) (eq. 5.7 of the paper), through two phases:
///
///
/// # Arguments
/// * `qtm`    — the Qtm being simulated; its transition function δ defines G
/// * `qstate` — register quantum state at the start of this step
///
/// # Returns
/// The register quantum state after one complete Qtm step
/// Returns an empty state immediately if `qstate` is empty
pub fn circuit_step(qtm: &Qtm, qstate: &RegQState) -> RegQState {
  let n = match qstate.0.keys().next() {
    Some(rs) => rs.n(),
    None => return RegQState::default(),
  };

  let mut state = qstate.clone();
  for centre in 0..n {
    state = apply_g_at(qtm, &state, centre);
  }

  state.f_all()
}

/// Runs the full circuit simulation for `steps` steps and returns the history
///
/// Sets up the tape loop, initialises the register state from `input`, and
/// calls [`circuit_step`] repeatedly, collecting the state after each step
///
////// # Arguments
/// * `qtm`   — the Qtm to simulate; must satisfy the Bernstein-Vazirani
///   conditions for the resulting circuit to be unitary
/// * `input` — tape symbols placed at positions 1..=|input|; must satisfy
///   |input| ≤ `steps` (the standard assumption from §5a)
/// * `steps` — number of simulation steps `t`; determines N and history
///   length
///
/// # Return value
/// A `Vec` of length `steps + 1` where:
/// * `history[0]`  is the initial state (before any step)
/// * `history[k]`  is the state after exactly `k` applications of
///   [`circuit_step`]
///
/// Keeping the full history allows the caller to animate the evolution or
/// compare every intermediate state against the direct simulation
pub fn circuit_run(qtm: &Qtm, input: &[usize], steps: usize) -> Vec<RegQState> {
  let n = tape_loop_size(steps);
  let init = RegState::from_input(input, n);
  let mut qstate = RegQState::from_reg_state(init);
  let mut history = vec![qstate.clone()]; // history[0] = initial state
  for _ in 0..steps {
    qstate = circuit_step(qtm, &qstate);
    history.push(qstate.clone());
  }
  history
}

/// Return the maximum |amplitude_direct − amplitude_circuit| across all
/// configurations reachable in `direct`.
pub fn verify_vs_direct(
  direct: &crate::qtm::QState,
  circuit: &RegQState,
  _t: usize,
) -> f64 {
  let n = match circuit.0.keys().next() {
    Some(rs) => rs.n(),
    None => return f64::NAN,
  };

  direct
    .0
    .iter()
    .map(|(cfg, &d_amp)| {
      // Build the register state that corresponds to this configuration.
      let mut rs = RegState::blank(n);
      for &(pos, sym) in &cfg.tape.0 {
        rs.0[pos.rem_euclid(n as i64) as usize].t = sym;
      }
      let head_idx = cfg.head.rem_euclid(n as i64) as usize;
      rs.0[head_idx].s = cfg.state as i32;

      let c_amp = circuit
        .0
        .get(&rs)
        .copied()
        .unwrap_or(crate::complex::C64::ZERO);
      (d_amp - c_amp).norm()
    })
    .fold(0.0_f64, f64::max)
}
