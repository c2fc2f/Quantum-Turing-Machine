# Quantum Turing Machine

A command-line tool written in Rust that simulates Quantum Turing Machines (QTMs) and converts them to equivalent quantum circuits using the method of Molina & Watrous (2019).

## Overview

The simulation runs in two parallel modes and cross-checks their results:

**Direct simulation** applies the global evolution operator U_δ (eq. 3.6) step by step to a sparse superposition of configurations. Each step branches each basis state into its transition targets, accumulating complex amplitudes that may interfere.

**Circuit simulation** implements the same evolution as a sequence of local quantum gates. Each step consists of two phases:
1. Gate G is applied at every tape-square position. G is a unitary acting on three consecutive cell-pairs derived from the QTM transition function via theorem 4.4 (localisation of causal unitary evolutions). The N applications per step commute on the valid-computation subspace im(Π) and can be parallelised into 3 depth-1 sublayers.
2. A global sign-flip F is applied to every state register, converting inactive shadow heads back to active heads.

The two simulations are verified to agree at machine-epsilon precision after every run.

**Reference:** Molina A, Watrous J. "Revisiting the simulation of quantum Turing machines by quantum circuits." *Proc. R. Soc. A* 475:20180767 (2019). https://doi.org/10.1098/rspa.2018.0767

## Requirements

- Rust toolchain (edition 2024, stable)
- Or: Nix with flakes enabled

## Installation

### From source

```
git clone https://github.com/c2fc2f/Quantum-Turing-Machine
cd Quantum-Turing-Machine
cargo build --release
```

The compiled binary will be at `target/release/qtm`.

### With Nix

A Nix flake is provided. Shell completions and man pages are installed automatically.

```
nix run github:c2fc2f/Quantum-Turing-Machine
# or
nix build
# or, to enter a development shell:
nix develop
```

## Usage

```
qtm
```

Running the binary without arguments executes all four built-in examples in sequence and prints a complexity summary.

### Shell completions

```
qtm completion bash  >> ~/.bashrc
qtm completion zsh   >> ~/.zshrc
qtm completion fish  >> ~/.config/fish/completions/qtm.fish
```

## Built-in examples

| Example | States | Symbols | Description |
|---|---|---|---|
| Shift-right | 1 | 2 | Deterministic head shift; validates infrastructure |
| Quantum walk | 2 | 2 | Hadamard-coin walk; genuinely quantum |
| Bit-flip | 1 | 2 | Classical symbol inversion |
| Flip-bounce | 2 | 2 | Flips and moves right in state 1, passes and moves left in state 2 |

Each example:
1. Checks the Bernstein-Vazirani unitarity conditions.
2. Runs both simulations.
3. Prints the full step-by-step history with amplitudes and head-position probabilities.
4. Reports the maximum amplitude error between the two simulations.

## Output format

Each cell `(Sᵢ, Tᵢ)` is displayed as:

```
[+q,s]   active head at this square, QTM state q, tape symbol s
[-q,s]   inactive (shadow) head, state q, tape symbol s
[ 0,s]   no head, tape symbol s
```

Head positions are reported in signed integer tape coordinates centred on 0, converted from the internal tape-loop index in Z_N.

## Complexity

For t steps on an input of length n ≤ t, choosing the tape loop size N as the smallest multiple of 3 ≥ 2t + 1:

| t | N | G-gates/step | Total G-gates |
|---|---|---|---|
| 1 | 3 | 3 | 3 |
| 4 | 9 | 9 | 36 |
| 16 | 33 | 33 | 528 |
| 64 | 129 | 129 | 8 256 |
| 128 | 258 | 258 | 33 024 |

- **Depth** O(t): 3 sublayers per step × t steps. Linear in t — the main improvement over Yao's original O(t²)-depth simulation.
- **Size** O(t²): N × t gates. Same asymptotic cost as Yao.

## Module structure

| Module | Contents |
|---|---|
| `complex` | 64-bit complex arithmetic; no external dependencies |
| `qtm` | QTM definition, Bernstein-Vazirani validity check, direct simulation |
| `register` | Cell, RegState, RegQState — register types for the circuit simulation |
| `gate_g` | Gate G implementation (equations 5.13–5.15) |
| `circuit` | Circuit simulation, tape-loop sizing, cross-verification |
| `cli` | Clap-based CLI definition and shell-completion subcommand |

## License

This project is licensed under the [MIT License](LICENSE).
