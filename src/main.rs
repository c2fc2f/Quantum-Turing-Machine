//! Binary to simulates Quantum Turing Machines (QTMs) and converts them to
//! equivalent quantum circuits

use std::{io::stdout, process::ExitCode};

use clap::{CommandFactory, Parser};
use clap_complete::generate;

use crate::cli::{Args, Command};

/// Binary to simulates Quantum Turing Machines (QTMs) and converts them to
/// equivalent quantum circuits
pub mod cli;

fn main() -> ExitCode {
  let args = Args::parse();

  match args.command {
    Some(Command::Completion(a)) => {
      let mut cmd = cli::Args::command();
      let name = cmd.get_name().to_string();
      generate(a.shell, &mut cmd, name, &mut stdout());
      ExitCode::SUCCESS
    }
    _ => ExitCode::SUCCESS,
  }
}
