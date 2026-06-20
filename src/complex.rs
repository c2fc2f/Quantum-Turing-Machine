//! Minimal 64-bit complex number arithmetic
//!
//! Provides [`C64`], a plain `(re, im): (f64, f64)` complex number with the
//! arithmetic operators and formatting needed by the QTM simulator. No
//! external crates are used; the implementation is intentionally small so
//! that every operation is easy to audit against the paper's formulae
//!
//! # Conventions
//! * Arithmetic follows the standard field rules over **\mathbb{C}**.
//! * "Near-zero" comparisons (e.g. amplitude pruning) are the caller's
//!   responsibility; this module only provides an exact [`C64::is_zero`]
//!   check
//! * [`fmt::Display`] and [`fmt::Debug`] produce the same output: the
//!   shortest f64 representation that round-trips, with zero components
//!   omitted

use std::{
  fmt,
  ops::{Add, AddAssign, Mul, Neg, Sub},
};

/// A complex number stored as two `f64` components
#[derive(Clone, Copy, PartialEq)]
pub struct C64 {
  /// Real part
  pub re: f64,
  /// Imaginary part
  pub im: f64,
}

impl C64 {
  /// The additive identity `0 + 0i`
  pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };

  /// The multiplicative identity `1 + 0i`
  pub const ONE: C64 = C64 { re: 1.0, im: 0.0 };

  /// The additive identity `\frac{1}{\sqrt{2}} + 0i`
  pub const INV_SQRT2: C64 = C64 {
    re: std::f64::consts::FRAC_1_SQRT_2,
    im: 0.0,
  };

  /// Constructs `re + im i`
  pub fn new(re: f64, im: f64) -> Self {
    C64 { re, im }
  }

  /// Constructs a purely real number `re + 0i`
  #[allow(dead_code)]
  pub fn real(re: f64) -> Self {
    C64 { re, im: 0.0 }
  }

  /// Returns the complex conjugate `re − im i`
  pub fn conj(self) -> Self {
    C64 {
      re: self.re,
      im: -self.im,
    }
  }

  /// Returns `{|z|}^2 = {re}^2 + {im}^2`
  ///
  /// Prefer this over [`norm`](Self::norm) whenever only the squared
  /// magnitude is needed, as it avoids the cost of a square-root.
  /// Amplitude pruning and probability calculations both use this form
  pub fn norm_sq(self) -> f64 {
    self.re * self.re + self.im * self.im
  }

  /// Returns `|z| = \sqrt{{re}^2 + {im}^2}`
  ///
  /// Use [`norm_sq`](Self::norm_sq) instead when comparing magnitudes,
  /// since squaring is monotone and cheaper
  pub fn norm(self) -> f64 {
    self.norm_sq().sqrt()
  }

  /// Returns `true` iff both components are exactly `0.0` at the bit level
  pub fn is_zero(self) -> bool {
    self == C64::ZERO
  }
}

impl Add for C64 {
  type Output = C64;
  fn add(self, r: C64) -> C64 {
    C64::new(self.re + r.re, self.im + r.im)
  }
}

impl AddAssign for C64 {
  fn add_assign(&mut self, r: C64) {
    self.re += r.re;
    self.im += r.im;
  }
}

impl Sub for C64 {
  type Output = C64;
  fn sub(self, r: C64) -> C64 {
    C64::new(self.re - r.re, self.im - r.im)
  }
}

impl Mul for C64 {
  type Output = C64;
  fn mul(self, r: C64) -> C64 {
    C64::new(
      self.re * r.re - self.im * r.im,
      self.re * r.im + self.im * r.re,
    )
  }
}

impl Neg for C64 {
  type Output = C64;
  fn neg(self) -> C64 {
    C64::new(-self.re, -self.im)
  }
}

impl fmt::Display for C64 {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match (self.re == 0.0, self.im == 0.0) {
      (true, true) => write!(f, "0"),
      (true, false) => write!(f, "{}i", self.im),
      (false, true) => write!(f, "{}", self.re),
      (false, false) if self.im >= 0.0 => write!(f, "{}+{}i", self.re, self.im),
      (false, false) => write!(f, "{}{}i", self.re, self.im),
    }
  }
}

impl fmt::Debug for C64 {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    fmt::Display::fmt(self, f)
  }
}
