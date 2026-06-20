//! Gate G — the local 3-cell unitary derived from the Qtm transition function.
//!
//! Reference: §5c-d of Molina & Watrous (2019), equations (5.13)–(5.15).
//!
//! All amplitude tests use **exact zero** (`== C64::ZERO`).
//! Tiny floating-point residuals are preserved at full IEEE-754 precision.

use crate::{complex::C64, qtm::Qtm, register::Cell};
use std::collections::HashMap;

/// Local 3-cell state: [left, centre, right].
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Local3(pub [Cell; 3]);

// ─────────────────────────────────────────────────────────────────────────────
// Apply gate G
// ─────────────────────────────────────────────────────────────────────────────

/// Apply gate G to a 3-cell local state; returns all (output, amplitude) branches.
///
/// Implements eqs. (5.13)–(5.15).  Only branches with **exactly** non-zero
/// amplitude are emitted; all others are naturally absent from the sum.
pub fn apply_g(qtm: &Qtm, local: Local3) -> Vec<(Local3, C64)> {
  let [c0, c1, c2] = local.0;

  let nz = [c0, c1, c2].iter().filter(|c| !c.is_empty()).count();

  // ── Trivial: 0 non-zero cells (empty neighbourhood)
  //            or ≥2 non-zero cells (multi-head, invalid) → identity ──────
  if nz != 1 {
    return vec![(local, C64::ONE)];
  }

  let (a1, a2, a3) = (c0.t, c1.t, c2.t);

  match (c0.s, c1.s, c2.s) {
    // Active head at centre → trivial (F leaves it alone, V cancels V*).
    (0, s, 0) if s > 0 => vec![(local, C64::ONE)],

    // Inactive head at left or right → trivial.
    (s, 0, 0) if s < 0 => vec![(local, C64::ONE)],
    (0, 0, s) if s < 0 => vec![(local, C64::ONE)],

    // ── Eq. (5.13): inactive head at centre ─────────────────────────────
    //
    // G |0,a₁⟩ |−p₂,a₂⟩ |0,a₃⟩  →
    //     Σ_{p₁,b₁} δ(p₁,b₁)[p₂,a₁,+1]  |p₁,b₁⟩ |0,a₂⟩ |0,a₃⟩
    //   + Σ_{p₃,b₃} δ(p₃,b₃)[p₂,a₃,−1]  |0,a₁⟩  |0,a₂⟩ |p₃,b₃⟩
    //
    // Meaning: the inactive head (state p₂) resulted from either
    //   (a) a head at the left cell moving right, or
    //   (b) a head at the right cell moving left.
    // G reconstructs the pre-step active head in each branch.
    (0, s1, 0) if s1 < 0 => {
      let p2 = c1.head_state();
      let mut out = Vec::new();

      // (a) head came from left — δ(p₁, b₁)[p₂, a₁, +1]
      for p1 in 1..=qtm.num_states {
        for b1 in 0..qtm.num_symbols {
          let amp = qtm.amp(p1, b1, p2, a1, 1);
          if !amp.is_zero() {
            out.push((
              Local3([
                Cell::active(p1, b1),
                Cell { s: 0, t: a2 },
                Cell { s: 0, t: a3 },
              ]),
              amp,
            ));
          }
        }
      }

      // (b) head came from right — δ(p₃, b₃)[p₂, a₃, −1]
      for p3 in 1..=qtm.num_states {
        for b3 in 0..qtm.num_symbols {
          let amp = qtm.amp(p3, b3, p2, a3, -1);
          if !amp.is_zero() {
            out.push((
              Local3([
                Cell { s: 0, t: a1 },
                Cell { s: 0, t: a2 },
                Cell::active(p3, b3),
              ]),
              amp,
            ));
          }
        }
      }
      out
    }

    // ── Eq. (5.14): active head at left ─────────────────────────────────
    //
    // G |p₁,a₁⟩ |0,a₂⟩ |0,a₃⟩  →
    //     Σ_{q₂,b₁}       δ(p₁,a₁)[q₂,b₁,+1]              |0,b₁⟩  |−q₂,a₂⟩ |0,a₃⟩
    //   + Σ_{q₁,b₁}  A[q₁,b₁]                              |q₁,b₁⟩ |0,a₂⟩  |0,a₃⟩
    //
    // where  A[q₁,b₁] = Σ_{r₀,c₁} δ(p₁,a₁)[r₀,c₁,−1] · δ(q₁,b₁)[r₀,c₁,−1]
    //   (product, NOT conjugate-product — derived from V*·Fᵢ·V via theorem 4.4)
    //
    // First sum:  head moves right → becomes inactive at centre.
    // Second sum: "bounce-back" through leftward V + rightward V* .
    //             Theorem 4.4 guarantees this is zero for valid Qtms;
    //             we compute it anyway and let the algebra decide.
    (s0, 0, 0) if s0 > 0 => {
      let p1 = c0.head_state();
      let mut out = Vec::new();

      // First sum
      for q2 in 1..=qtm.num_states {
        for b1 in 0..qtm.num_symbols {
          let amp = qtm.amp(p1, a1, q2, b1, 1);
          if !amp.is_zero() {
            out.push((
              Local3([
                Cell { s: 0, t: b1 },
                Cell::inactive(q2, a2),
                Cell { s: 0, t: a3 },
              ]),
              amp,
            ));
          }
        }
      }

      // Second sum: accumulate A[q₁,b₁] = Σ_{r₀,c₁} δ(p₁,a₁)[r₀,c₁,−1]·δ(q₁,b₁)[r₀,c₁,−1]
      let mut coeff: HashMap<(usize, usize), C64> = HashMap::new();
      for r0 in 1..=qtm.num_states {
        for c1s in 0..qtm.num_symbols {
          let ap1 = qtm.amp(p1, a1, r0, c1s, -1);
          if ap1.is_zero() {
            continue;
          }
          for q1 in 1..=qtm.num_states {
            for b1 in 0..qtm.num_symbols {
              let aq1 = qtm.amp(q1, b1, r0, c1s, -1);
              if !aq1.is_zero() {
                *coeff.entry((q1, b1)).or_insert(C64::ZERO) += ap1 * aq1;
              }
            }
          }
        }
      }
      for ((q1, b1), amp) in coeff {
        if !amp.is_zero() {
          out.push((
            Local3([
              Cell::active(q1, b1),
              Cell { s: 0, t: a2 },
              Cell { s: 0, t: a3 },
            ]),
            amp,
          ));
        }
      }
      out
    }

    // ── Eq. (5.15): active head at right ────────────────────────────────
    //
    // G |0,a₁⟩ |0,a₂⟩ |p₃,a₃⟩  →
    //     Σ_{q₂,b₃}       δ(p₃,a₃)[q₂,b₃,−1]              |0,a₁⟩ |−q₂,a₂⟩ |0,b₃⟩
    //   + Σ_{q₃,b₃}  B[q₃,b₃]                              |0,a₁⟩ |0,a₂⟩  |q₃,b₃⟩
    //
    // where  B[q₃,b₃] = Σ_{r₄,c₃} δ(p₃,a₃)[r₄,c₃,+1] · δ(q₃,b₃)[r₄,c₃,+1]
    //
    // Mirror image of eq. (5.14) for the rightward case.
    (0, 0, s2) if s2 > 0 => {
      let p3 = c2.head_state();
      let mut out = Vec::new();

      // First sum
      for q2 in 1..=qtm.num_states {
        for b3 in 0..qtm.num_symbols {
          let amp = qtm.amp(p3, a3, q2, b3, -1);
          if !amp.is_zero() {
            out.push((
              Local3([
                Cell { s: 0, t: a1 },
                Cell::inactive(q2, a2),
                Cell { s: 0, t: b3 },
              ]),
              amp,
            ));
          }
        }
      }

      // Second sum: B[q₃,b₃] = Σ_{r₄,c₃} δ(p₃,a₃)[r₄,c₃,+1]·δ(q₃,b₃)[r₄,c₃,+1]
      let mut coeff: HashMap<(usize, usize), C64> = HashMap::new();
      for r4 in 1..=qtm.num_states {
        for c3s in 0..qtm.num_symbols {
          let ap3 = qtm.amp(p3, a3, r4, c3s, 1);
          if ap3.is_zero() {
            continue;
          }
          for q3 in 1..=qtm.num_states {
            for b3 in 0..qtm.num_symbols {
              let aq3 = qtm.amp(q3, b3, r4, c3s, 1);
              if !aq3.is_zero() {
                *coeff.entry((q3, b3)).or_insert(C64::ZERO) += ap3 * aq3;
              }
            }
          }
        }
      }
      for ((q3, b3), amp) in coeff {
        if !amp.is_zero() {
          out.push((
            Local3([
              Cell { s: 0, t: a1 },
              Cell { s: 0, t: a2 },
              Cell::active(q3, b3),
            ]),
            amp,
          ));
        }
      }
      out
    }

    _ => vec![(local, C64::ONE)],
  }
}
