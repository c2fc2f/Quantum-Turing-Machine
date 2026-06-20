//! Qtm model (§3, Molina & Watrous 2019) and direct simulation.
//!
//! All amplitude filtering uses **exact zero** only — no epsilon rounding.
//! The Bernstein-Vazirani validity check uses a small tolerance (1e-8)
//! because it compares floating-point sums against the mathematical target
//! values 0 and 1; that is a correctness assertion, not an amplitude-pruning
//! step.

use crate::complex::C64;
use std::{collections::HashMap, fmt};

/// One output branch of δ(p, a): the amplitude to transition to state `new_state`,
/// write `write_sym`, and move in `direction` ∈ {−1, +1}.
#[derive(Clone, Debug)]
pub struct TEntry {
  /// Destination Qtm state q ∈ {1, …, m}.
  pub new_state: usize,
  /// Symbol written to the current tape square.
  pub write_sym: usize,
  /// Head movement: −1 (left) or +1 (right).
  pub direction: i32,
  /// Transition amplitude δ(p,a)[q,b,D].
  pub amplitude: C64,
}

/// A Quantum Turing Machine.
///
/// State set Q = {1, …, `num_states`}; tape alphabet Γ = {0, …, `num_symbols`−1}
/// where 0 is the blank.  Transition amplitudes are stored in `delta` as
/// `(state, read_symbol) → [TEntry]`; see [`add_transition`](Self::add_transition).
pub struct Qtm {
  /// Number of states m; valid states are 1..=num_states.
  pub num_states: usize,
  /// Size of the tape alphabet k; symbols are 0..num_symbols (0 = blank).
  pub num_symbols: usize,
  /// Transition function: (from_state, read_symbol) → output branches.
  pub delta: HashMap<(usize, usize), Vec<TEntry>>,
}

impl Qtm {
  /// Creates an empty Qtm with the given state and symbol counts.
  ///
  /// Panics if `num_states < 1` or `num_symbols < 2` (need at least blank + one symbol).
  pub fn new(num_states: usize, num_symbols: usize) -> Self {
    assert!(num_states >= 1);
    assert!(num_symbols >= 2);
    Qtm {
      num_states,
      num_symbols,
      delta: HashMap::new(),
    }
  }

  /// Registers δ(p, a) as a list of `(new_state, write_sym, direction, amplitude)` tuples.
  ///
  /// Entries with amplitude exactly 0.0 are silently dropped at insertion time.
  pub fn add_transition(
    &mut self,
    p: usize,
    a: usize,
    entries: &[(usize, usize, i32, C64)],
  ) {
    let te: Vec<TEntry> = entries
      .iter()
      .filter(|&&(_, _, _, amp)| !amp.is_zero())
      .map(|&(q, b, d, amp)| TEntry {
        new_state: q,
        write_sym: b,
        direction: d,
        amplitude: amp,
      })
      .collect();
    self.delta.insert((p, a), te);
  }

  /// Returns δ(p, a)[q, b, D], or exactly 0.0 if the branch is not defined.
  pub fn amp(&self, p: usize, a: usize, q: usize, b: usize, d: i32) -> C64 {
    self
      .delta
      .get(&(p, a))
      .and_then(|es| {
        es.iter()
          .find(|e| e.new_state == q && e.write_sym == b && e.direction == d)
      })
      .map_or(C64::ZERO, |e| e.amplitude)
  }

  /// Checks both Bernstein-Vazirani conditions (tolerance 1e-8 on inner-product sums).
  ///
  /// * **Condition 1**: {δ(p,a) : p ∈ Q, a ∈ Γ} is orthonormal.
  /// * **Condition 2**: Σ_q δ(p₀,a₀)[q,b₀,+1] · conj(δ(p₁,a₁)[q,b₁,−1]) = 0 for all inputs.
  ///
  /// Returns `Ok(())` if both hold, or `Err` with a description of the first failure found.
  pub fn is_valid(&self) -> Result<(), String> {
    let eps = 1e-8;
    for p0 in 1..=self.num_states {
      for a0 in 0..self.num_symbols {
        for p1 in 1..=self.num_states {
          for a1 in 0..self.num_symbols {
            let mut dot = C64::ZERO;
            for q in 1..=self.num_states {
              for b in 0..self.num_symbols {
                for &d in &[-1i32, 1] {
                  dot += self.amp(p0, a0, q, b, d)
                    * self.amp(p1, a1, q, b, d).conj();
                }
              }
            }
            let want = if p0 == p1 && a0 == a1 { 1.0 } else { 0.0 };
            if (dot.re - want).abs() > eps || dot.im.abs() > eps {
              return Err(format!(
                "Cond.1 failed: <δ({p0},{a0})|δ({p1},{a1})> = {dot} (want {want})"
              ));
            }
          }
        }
      }
    }
    for p0 in 1..=self.num_states {
      for a0 in 0..self.num_symbols {
        for b0 in 0..self.num_symbols {
          for p1 in 1..=self.num_states {
            for a1 in 0..self.num_symbols {
              for b1 in 0..self.num_symbols {
                let mut s = C64::ZERO;
                for q in 1..=self.num_states {
                  s += self.amp(p0, a0, q, b0, 1)
                    * self.amp(p1, a1, q, b1, -1).conj();
                }
                if s.norm_sq() > eps * eps {
                  return Err(format!(
                    "Cond.2 failed ({p0},{a0},{b0}),({p1},{a1},{b1}): {s}"
                  ));
                }
              }
            }
          }
        }
      }
    }
    Ok(())
  }
}

/// Sparse integer-indexed tape; only non-blank (non-zero) positions are stored, sorted by key.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct Tape(
  /// Sorted list of `(position, symbol)` pairs; blank cells (symbol 0) are absent.
  pub Vec<(i64, usize)>,
);

impl Tape {
  /// Creates an empty tape (all squares blank).
  pub fn new() -> Self {
    Tape(Vec::new())
  }

  /// Returns the symbol at `pos`, or 0 (blank) if absent.
  pub fn get(&self, pos: i64) -> usize {
    self
      .0
      .binary_search_by_key(&pos, |&(p, _)| p)
      .map(|i| self.0[i].1)
      .unwrap_or(0)
  }

  /// Returns a new tape with position `pos` set to `sym`.
  /// Writing 0 (blank) removes the entry.
  pub fn set(&self, pos: i64, sym: usize) -> Self {
    let mut v = self.0.clone();
    match v.binary_search_by_key(&pos, |&(p, _)| p) {
      Ok(i) => {
        if sym == 0 {
          v.remove(i);
        } else {
          v[i] = (pos, sym);
        }
      }
      Err(i) => {
        if sym != 0 {
          v.insert(i, (pos, sym));
        }
      }
    }
    Tape(v)
  }

  /// Builds a tape from an input slice: `input[i]` lands at position `i + 1`.
  /// Blank symbols (0) are not stored.
  pub fn from_input(input: &[usize]) -> Self {
    let mut t = Tape::new();
    for (i, &s) in input.iter().enumerate() {
      if s != 0 {
        t = t.set(i as i64 + 1, s);
      }
    }
    t
  }
}

impl fmt::Display for Tape {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{{")?;
    for (i, &(pos, sym)) in self.0.iter().enumerate() {
      if i > 0 {
        write!(f, ", ")?;
      }
      write!(f, "{pos}↦{sym}")?;
    }
    write!(f, "}}")
  }
}

/// A classical Qtm configuration: (state, head position, tape contents).
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Config {
  /// Current Qtm state ∈ {1, …, m}.
  pub state: usize,
  /// Integer tape-head position.
  pub head: i64,
  /// Tape contents (only non-blank cells stored).
  pub tape: Tape,
}

impl Config {
  /// Initial configuration for `input`: state 1, head at 0, input symbols at positions 1..n.
  pub fn initial(input: &[usize]) -> Self {
    Config {
      state: 1,
      head: 0,
      tape: Tape::from_input(input),
    }
  }
}

impl fmt::Display for Config {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(
      f,
      "(q={}, h={:+}, tape={})",
      self.state, self.head, self.tape
    )
  }
}

/// A quantum superposition of Qtm configurations: Σ amplitude[config] |config⟩.
#[derive(Clone, Default)]
pub struct QState(
  /// Sparse map from basis configurations to complex amplitudes.
  pub HashMap<Config, C64>,
);

impl QState {
  /// Constructs the pure state `|c⟩` with amplitude 1.
  pub fn from_config(c: Config) -> Self {
    let mut s = QState::default();
    s.0.insert(c, C64::ONE);
    s
  }

  /// Adds `a` to the amplitude of configuration `c`.
  pub fn add_amp(&mut self, c: Config, a: C64) {
    *self.0.entry(c).or_insert(C64::ZERO) += a;
  }

  /// Removes entries with amplitude exactly 0.0 — no epsilon.
  pub fn clean(&mut self) {
    self.0.retain(|_, a| !a.is_zero());
  }

  /// Returns ‖ψ‖² = Σ |amplitude|²; should equal 1.0 for a valid Qtm.
  pub fn total_prob(&self) -> f64 {
    self.0.values().map(|a| a.norm_sq()).sum()
  }
}

impl fmt::Display for QState {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    let mut cfgs: Vec<_> = self.0.iter().collect();
    cfgs.sort_by_key(|(c, _)| (c.head, c.state));
    for (cfg, amp) in &cfgs {
      writeln!(f, "      {amp:>30}  ×  {cfg}")?;
    }
    Ok(())
  }
}

/// Applies U_δ once: |p, i, T⟩ → Σ_{q,b,D} δ(p,T(i))[q,b,D] |q, i+D, T_{i←b}⟩.
///
/// Amplitudes accumulated from different source configurations may interfere.
/// The result is cleaned of exact-zero entries before being returned.
pub fn direct_step(qtm: &Qtm, state: &QState) -> QState {
  let mut next = QState::default();
  for (cfg, &amp) in &state.0 {
    let sym = cfg.tape.get(cfg.head);
    if let Some(entries) = qtm.delta.get(&(cfg.state, sym)) {
      for e in entries {
        // Exact-zero amplitudes were already filtered at add_transition time.
        next.add_amp(
          Config {
            state: e.new_state,
            head: cfg.head + e.direction as i64,
            tape: cfg.tape.set(cfg.head, e.write_sym),
          },
          amp * e.amplitude,
        );
      }
    }
  }
  next.clean();
  next
}

/// Runs the direct simulation for `steps` steps, returning the full history.
///
/// `history[0]` is the initial state; `history[k]` is the state after `k` applications
/// of [`direct_step`].  The vector is pre-allocated to `steps + 1` entries.
pub fn direct_run(qtm: &Qtm, input: &[usize], steps: usize) -> Vec<QState> {
  let mut history = Vec::with_capacity(steps + 1);
  let mut state = QState::from_config(Config::initial(input));
  history.push(state.clone());
  for _ in 0..steps {
    state = direct_step(qtm, &state);
    history.push(state.clone());
  }
  history
}
