//! Register types for the circuit simulation (§5, Molina & Watrous 2019)
//!
//! Three types form a hierarchy:
//! * [`Cell`] — one tape square (Sᵢ, Tᵢ)
//! * [`RegState`] — a full N-cell classical configuration
//! * [`RegQState`] — a quantum superposition of `RegState`s

use crate::complex::C64;
use std::{collections::HashMap, fmt};

/// One tape-square register pair (Sᵢ, Tᵢ)
///
/// `s` encodes head presence and QTM state:
/// * `s > 0` — active head, QTM state `s`
/// * `s < 0` — inactive (shadow) head, state `|s|`
/// * `s = 0` — no head at this square
///
/// `t` is the tape symbol ∈ {0, …, k−1}, where 0 is the blank
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub struct Cell {
  /// Head indicator and QTM state; see type-level docs for the sign
  /// convention
  pub s: i32,
  /// Tape symbol ∈ {0, …, k−1}
  pub t: usize,
}

impl Cell {
  /// A blank cell with no head: `s = 0`, `t = 0`
  pub fn blank() -> Self {
    Cell { s: 0, t: 0 }
  }

  /// A cell carrying an **active** head in QTM state `state`, with symbol
  /// `sym`
  pub fn active(state: usize, sym: usize) -> Self {
    Cell {
      s: state as i32,
      t: sym,
    }
  }

  /// A cell carrying an **inactive** (shadow) head in state `state`, with
  /// symbol `sym`
  ///
  /// Inactive heads are created by gate G and converted back to active by
  /// the global F flip at the end of each simulation step
  pub fn inactive(state: usize, sym: usize) -> Self {
    Cell {
      s: -(state as i32),
      t: sym,
    }
  }

  /// F operator (eq. 5.6): flips the sign of `s`, leaving `t` unchanged
  ///
  /// Maps active ↔ inactive and leaves empty cells fixed (`0 ↦ 0`)
  pub fn f(self) -> Self {
    Cell {
      s: -self.s,
      t: self.t,
    }
  }

  /// Returns `true` if this cell holds an active head (`s > 0`)
  pub fn is_active(self) -> bool {
    self.s > 0
  }

  /// Returns `true` if no head is present at this cell (`s = 0`)
  pub fn is_empty(self) -> bool {
    self.s == 0
  }

  /// Returns the head state `|s|`, regardless of active/inactive sign
  ///
  /// Returns 0 for empty cells; callers should guard with
  /// [`is_active`](Self::is_active) or [`is_empty`](Self::is_empty) when the
  /// sign matters
  pub fn head_state(self) -> usize {
    self.s.unsigned_abs() as usize
  }
}

impl fmt::Display for Cell {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self.s.cmp(&0) {
      std::cmp::Ordering::Greater => write!(f, "[+{},{}]", self.s, self.t),
      std::cmp::Ordering::Less => write!(f, "[-{},{}]", -self.s, self.t),
      std::cmp::Ordering::Equal => write!(f, "[ 0,{}]", self.t),
    }
  }
}

/// A classical register configuration: N cell-pairs indexed by Z_N
///
/// The inner `Vec<Cell>` has length N throughout a simulation run; index `i`
/// corresponds to tape-loop position `i ∈ {0, …, N−1}`
#[derive(Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub struct RegState(pub Vec<Cell>);

impl RegState {
  /// Returns an all-blank N-cell state (no head, all symbols 0)
  pub fn blank(n: usize) -> Self {
    RegState(vec![Cell::blank(); n])
  }

  /// Returns the tape-loop length N
  pub fn n(&self) -> usize {
    self.0.len()
  }

  /// Builds the initial register state for a given QTM input string
  ///
  /// * Cell 0: active head in state 1 (`s = 1`, `t = 0`)
  /// * Cells 1..=|input|: `t = input[i−1]` (symbols at tape positions 1..n)
  /// * All other cells: blank
  pub fn from_input(input: &[usize], n: usize) -> Self {
    let mut rs = RegState::blank(n);
    rs.0[0].s = 1;
    for (i, &sym) in input.iter().enumerate() {
      rs.0[(i + 1) % n].t = sym;
    }
    rs
  }

  /// Applies F to every cell simultaneously (global F₀⋯F_{N-1}, eq. 5.7)
  pub fn f_all(&self) -> Self {
    RegState(self.0.iter().map(|c| c.f()).collect())
  }

  /// Returns `(index, state)` of the unique active head, or `None` if absent
  pub fn active_head(&self) -> Option<(usize, usize)> {
    self.0.iter().enumerate().find_map(|(i, c)| {
      if c.is_active() {
        Some((i, c.head_state()))
      } else {
        None
      }
    })
  }
}

impl fmt::Display for RegState {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    for (i, c) in self.0.iter().enumerate() {
      if i > 0 {
        write!(f, " ")?;
      }
      write!(f, "{c}")?;
    }
    Ok(())
  }
}

/// A quantum superposition of register states: Σ amplitude\[rs\] |rs⟩.
///
/// The map is sparse — only basis states with non-zero amplitude are stored.
/// Use [`add_amp`](Self::add_amp) to accumulate amplitudes and
/// [`clean`](Self::clean) to prune exact zeros after a gate application
#[derive(Clone, Default)]
pub struct RegQState(pub HashMap<RegState, C64>);

impl RegQState {
  /// Constructs the pure state `|rs⟩` with amplitude 1
  pub fn from_reg_state(rs: RegState) -> Self {
    let mut s = RegQState::default();
    s.0.insert(rs, C64::ONE);
    s
  }

  /// Adds `a` to the amplitude of `rs`, inserting the entry if absent
  pub fn add_amp(&mut self, rs: RegState, a: C64) {
    *self.0.entry(rs).or_insert(C64::ZERO) += a;
  }

  /// Removes entries whose amplitude is exactly `C64::ZERO` (bit-level)
  ///
  /// Call after each gate application to keep the map compact.
  /// This is an exact check; near-zero entries from floating-point
  /// rounding are not removed
  pub fn clean(&mut self) {
    self.0.retain(|_, a| !a.is_zero());
  }

  /// Returns the total probability ‖ψ‖² = Σ |amplitude|²
  ///
  /// Should equal 1.0 throughout a valid simulation; deviations indicate
  /// a non-unitary QTM or accumulated floating-point error
  pub fn total_prob(&self) -> f64 {
    self.0.values().map(|a| a.norm_sq()).sum()
  }

  /// Applies F to every S register of every basis state, returning a new
  /// superposition
  ///
  /// Implements the global F₀⋯F_{N-1} phase of one circuit step (§5b)
  pub fn f_all(&self) -> Self {
    let mut next = RegQState::default();
    for (rs, &amp) in &self.0 {
      next.add_amp(rs.f_all(), amp);
    }
    next
  }

  /// Returns the marginal probability distribution over active-head
  /// positions
  ///
  /// Iterates all basis states, finds the active head in each, and sums
  /// `|amplitude|²` by `(tape_index, QTM_state)` pair.  The result is
  /// sorted by `(index, state)` for stable output
  ///
  /// Basis states with no active head (e.g. mid-step intermediate states)
  /// are silently skipped
  pub fn head_distribution(&self) -> Vec<(usize, usize, f64)> {
    let mut map: HashMap<(usize, usize), f64> = HashMap::new();
    for (rs, amp) in &self.0 {
      if let Some((idx, st)) = rs.active_head() {
        *map.entry((idx, st)).or_insert(0.0) += amp.norm_sq();
      }
    }
    let mut v: Vec<_> = map.into_iter().map(|((i, s), p)| (i, s, p)).collect();
    v.sort_by_key(|&(i, s, _)| (i, s));
    v
  }
}

impl fmt::Display for RegQState {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    let mut states: Vec<_> = self.0.iter().collect();
    states.sort_by_key(|(rs, _)| (*rs).clone());
    for (rs, amp) in &states {
      writeln!(f, "      {amp:>30}  ×  [{rs}]")?;
    }
    Ok(())
  }
}
