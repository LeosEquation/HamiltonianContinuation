# HamiltonianContinuation.jl

[![Julia](https://img.shields.io/badge/Julia-1.12%2B-9558B2.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**HamiltonianContinuation.jl** is a Julia package for the numerical continuation and bifurcation analysis of equilibria and periodic orbits in Hamiltonian dynamical systems.

The package is designed for computational studies of nonlinear dynamical systems, with an emphasis on high-performance numerical continuation, stability analysis, and the detection and continuation of bifurcations.

---

## Features

HamiltonianContinuation provides tools for:

- Numerical continuation of equilibrium branches.
- Numerical continuation of periodic-orbit branches.
- Detection and continuation of equilibrium bifurcations.
- Detection and continuation of limit points.
- Detection and continuation of branch points.
- Detection and continuation of flip (period-doubling) bifurcations.
- Stability analysis through eigenvalues and Floquet multipliers.
- Multiple-shooting methods for periodic orbits.
- Taylor-series-based numerical integration.
- Poincaré sections and maps.
- Memory-efficient repeated integrations for periodic-orbit continuation.

The package is built around Julia's type system and is intended to support large continuation calculations where numerical objects are repeatedly updated in memory.

---

## Installation

From the Julia package manager:

```julia
using Pkg
Pkg.add("HamiltonianContinuation")
```

Or from the Julia REPL:

```text
] add HamiltonianContinuation
```

For the development version:

```julia
] add https://github.com/LeosEquation/HamiltonianContinuation.jl
```

---

## Basic usage

Load the package with:

```julia
using HamiltonianContinuation
```

The main continuation routines are:

```julia
EquilibriaContinuation(...)
OrbitContinuation(...)
LPContinuation(...)
BPContinuation(...)
FPContinuation(...)
```

where:

- `EquilibriaContinuation` performs continuation of equilibrium branches.
- `OrbitContinuation` performs continuation of periodic orbits.
- `LPContinuation` performs continuation from limit points.
- `BPContinuation` performs continuation from branch points.
- `FPContinuation` performs continuation from flip (period-doubling) points.

The continuation routines return typed objects containing the computed branch and the associated bifurcation information.

---

## Main data structures

HamiltonianContinuation separates the numerical continuation algorithms from the objects used to store their results.

### Equilibrium branches

```julia
HamiltonianEquilibriumBranch
```

stores the data associated with a continued equilibrium branch, including equilibrium coordinates, energy, stability information, and detected bifurcation points.

### Periodic-orbit branches

```julia
HamiltonianPeriodicBranch
```

stores the data associated with a continued family of periodic orbits, including orbit information, period, and Floquet stability data.

### Limit points

```julia
HamiltonianLimitPoint
HamiltonianLimitPointBranch
```

represent detected limit points and the branches continued from them.

### Branch points

```julia
HamiltonianBranchPoint
HamiltonianBranchPointBranch
```

represent detected branch points and the branches continued from them.

### Flip points

```julia
HamiltonianFlipPoint
HamiltonianFlipPointBranch
```

represent detected flip (period-doubling) bifurcations of periodic orbits and the branches continued from them.

---

## Numerical methods

### Equilibrium continuation

Equilibrium branches are computed using numerical continuation combined with Newton correction.

The stability of equilibria is determined from the eigenvalues of the linearized dynamical system. Singular points detected along a branch can be stored and subsequently used to initiate dedicated bifurcation continuation procedures.

### Periodic-orbit continuation

Periodic orbits are computed using a multiple-shooting formulation. The resulting nonlinear system is corrected using Newton's method.

The stability of periodic orbits is analyzed through the monodromy matrix and its Floquet multipliers.

For Hamiltonian systems, the symplectic structure imposes constraints on the Floquet spectrum. These properties can be exploited when analyzing the stability and bifurcations of periodic solutions.

### Taylor-series integration

HamiltonianContinuation uses Taylor-series-based numerical integration for the repeated integrations required by continuation and multiple-shooting algorithms.

Periodic-orbit continuation can require thousands of integrations in a single computation. For calculations where only the final state is required, specialized integration routines avoid storing intermediate integration steps and return the final result directly, substantially reducing memory usage and computational overhead.

---

## Poincaré sections

HamiltonianContinuation provides functionality for constructing and evaluating Poincaré sections and maps.

These tools can be used to investigate:

- periodic orbits,
- invariant structures,
- resonances,
- stability,
- bifurcations,
- and the local phase-space structure surrounding periodic solutions.

---

## Architecture

The package is organized according to the main dynamical objects and continuation procedures.

```text
HamiltonianContinuation/
├── Project.toml
├── Manifest.toml
├── src/
│   ├── HamiltonianContinuation.jl
│   │
│   ├── objects/
│   │   └── ...
│   │
│   ├── equilibria/
│   │   ├── equilibria.jl
│   │   ├── eqbranchswitching.jl
│   │   ├── eqcontinuation.jl
│   │   ├── eqfinding.jl
│   │   └── eqsystem.jl
│   │
│   ├── limitpoint/
│   │   └── ...
│   │
│   ├── branchpoint/
│   │   └── ...
│   │
│   ├── orbits/
│   │   └── ...
│   │
│   └── flippoint/
│       └── ...
│
└── test/
    └── ...
```

The source tree is divided into separate components for equilibria, periodic orbits, and bifurcation-specific continuation algorithms.

---

## Dependencies

HamiltonianContinuation relies on the following Julia packages:

- [TaylorIntegration.jl](https://github.com/PerezHz/TaylorIntegration.jl) — Taylor-series-based numerical integration.
- [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) — Taylor polynomial arithmetic.
- [PeriodicSchurDecompositions.jl](https://github.com/andreasnoack/PeriodicSchurDecompositions.jl) — periodic Schur decomposition and Floquet analysis.
- `LinearAlgebra` — linear algebra functionality from the Julia standard library.
- `Printf` — formatted output from the Julia standard library.

A customized development version of `TaylorIntegration.jl` is currently used for specialized high-performance integrations required by the periodic-orbit continuation algorithms.

The customized version is maintained in a fork of the original repository while these modifications are being developed and evaluated.

---

## Development

Clone the repository:

```bash
git clone https://github.com/LeosEquation/HamiltonianContinuation.jl.git
cd HamiltonianContinuation.jl
```

Activate and instantiate the Julia environment:

```julia
using Pkg

Pkg.activate(".")
Pkg.instantiate()
```

To run the test suite:

```julia
Pkg.test()
```

or from the package manager:

```text
] test
```

---

## Development of TaylorIntegration

HamiltonianContinuation currently uses a customized development version of `TaylorIntegration.jl`.

The modifications primarily target repeated integrations required by periodic-orbit continuation, where storing intermediate integration steps can represent a significant computational and memory overhead.

The customized implementation provides specialized functionality for integrations where only the final state is required, together with additional functionality needed by the continuation algorithms.

These modifications are maintained separately from HamiltonianContinuation to preserve the separation between the numerical integration infrastructure and the continuation algorithms themselves.

---

## Scientific scope

HamiltonianContinuation is intended for computational studies of nonlinear Hamiltonian systems, particularly problems where continuation requires a large number of high-accuracy numerical integrations.

Potential applications include:

- Hamiltonian dynamical systems.
- Nonlinear dynamics.
- Numerical bifurcation analysis.
- Equilibrium continuation.
- Periodic-orbit continuation.
- Floquet stability analysis.
- Poincaré maps and sections.
- Resonances and invariant structures.
- Phase-space exploration.

The package was initially developed in the context of low-dimensional Hamiltonian systems and is being designed with extensibility toward broader continuation problems in mind.

---

## Current status

HamiltonianContinuation is currently under active development.

The core numerical continuation algorithms are being organized into a reusable Julia package, while the public API, documentation, tests, and examples are progressively being expanded.

The current implementation is primarily focused on Hamiltonian systems with low-dimensional reduced phase spaces, with particular emphasis on high-accuracy and high-performance repeated integrations.

---

## Citation

If you use HamiltonianContinuation in academic work, please cite the associated research work:

```bibtex
@software{HamiltonianContinuation,
  author  = {Mayorga López, Leonel},
  title   = {HamiltonianContinuation.jl},
  year    = {2026},
  url     = {https://github.com/LeosEquation/HamiltonianContinuation.jl}
}
```

A formal citation entry will be provided once the associated research work is published.

---

## License

HamiltonianContinuation.jl is released under the MIT License.

See [LICENSE](LICENSE) for details.