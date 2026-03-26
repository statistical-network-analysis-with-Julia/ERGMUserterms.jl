# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ERGMUserterms.jl is a Julia port of the R `ergm.userterms` package (StatNet collection) that provides templates, utilities, and validation tools for developing custom ERGM (Exponential Random Graph Model) terms. It integrates with ERGM.jl and Network.jl for model estimation.

## Development Commands

- **Run tests:** `julia --project -e 'using Pkg; Pkg.test()'`
- **Load package in REPL:** `julia --project -e 'using ERGMUserterms'`
- **Build docs:** `julia --project=docs docs/make.jl`

## Architecture

The package is a single-module design in `src/ERGMUserterms.jl` with no sub-files. It contains:

- **Abstract type:** `AbstractUserTerm <: AbstractERGMTerm` -- base type for all user-defined terms.
- **Term interface:** Every term must implement three methods: `name(term) -> String`, `compute(term, net) -> Float64` (full network statistic), and `change_stat(term, net, i, j) -> Float64` (change when toggling edge (i,j)). The invariant is `change_stat(term, net, i, j) == compute(net_toggled) - compute(net)`.
- **`@ergm_term` macro:** Convenience macro for defining terms with automatic validation hooks.
- **Validation/testing:** `validate_term()`, `change_stat_check()`, `consistency_check()`, `test_term()` -- check correctness and consistency of term implementations.
- **Benchmarking:** `benchmark_term()` profiles `compute()` vs `change_stat()` performance.
- **Example terms:** `ExampleTerm`, `TemplateTerm`, `WeightedEdges`, `DyadCovTerm`, `InteractionTerm` -- shipped as reference implementations.
- **Doc helpers:** `term_signature()`, `term_documentation()` for auto-generating term docs.

## Key Dependencies

- **ERGM.jl** -- provides `AbstractERGMTerm`, core ERGM functionality
- **Network.jl** -- network data structure (`Network`, `has_edge`, `add_edge!`, `rem_edge!`, `nv`, `ne`, etc.)
- **Graphs.jl** -- graph algorithms and iteration (`edges`, `vertices`, `src`, `dst`, `outneighbors`, `inneighbors`)
- Requires Julia >= 1.9

## Conventions

- Term structs should be **immutable** (`struct`, not `mutable struct`) and subtype `AbstractUserTerm`.
- `change_stat()` should be O(degree), not O(edges), for performance.
- Term names should be descriptive and include parameters (e.g., `"template.$(t.param)"`).
- Validation uses a tolerance of `1e-10` for floating-point consistency checks.
- All exported symbols are declared at the top of the module; there are no separate include files.
