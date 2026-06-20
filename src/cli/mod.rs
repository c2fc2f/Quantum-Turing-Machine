//! Command-line interface definition for the Quantum-Turing-Machine

pub mod completion;

use clap::{Parser, Subcommand};

/// Simulates Quantum Turing Machines (QTMs) and converts them to equivalent
/// quantum circuits using the method of Molina & Watrous (2019).
///
/// Reference:
///   Molina A, Watrous J. "Revisiting the simulation of quantum Turing
///   machines by quantum circuits." Proc. R. Soc. A 475:20180767 (2019).
///   https://doi.org/10.1098/rspa.2018.0767
#[derive(Parser, Debug)]
#[command(version, about, long_about)]
pub struct Args {
  /// The specific operation to perform with the binary
  #[command(subcommand)]
  pub command: Option<Command>,
}

/// List of available subcommands in the binary
#[derive(Subcommand, Debug)]
#[non_exhaustive]
pub enum Command {
  /// Print shell completions and exit
  #[command(hide = true)]
  Completion(completion::Args),
}
