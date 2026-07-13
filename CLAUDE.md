# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ERGMUserterms.jl is a Julia port of the R `ergm.userterms` package (StatNet collection) that provides templates, utilities, and validation tools for developing custom ERGM (Exponential Random Graph Model) terms. It integrates with ERGM.jl and Networks.jl for model estimation.

## Development Commands

- **Run tests:** `julia --project -e 'using Pkg; Pkg.test()'`
- **Load package in REPL:** `julia --project -e 'using ERGMUserterms'`
- **Build docs:** `julia --project=docs docs/make.jl`

## Architecture

The package is a single-module design in `src/ERGMUserterms.jl` with no sub-files. It contains:

- **Abstract type:** `AbstractUserTerm <: AbstractERGMTerm` -- base type for all user-defined terms.
- **Term interface:** Every term must implement three methods: `name(term) -> String`, `compute(term, net) -> Float64` (full network statistic), and `change_stat(term, net, i, j) -> Float64` — the **add-direction** change statistic `g(y⁺ij) − g(y⁻ij)` (statistic with edge (i,j) present minus absent), independent of the dyad's current state. The toggle-direction idiom `has_edge ? -Δ : Δ` is wrong and is rejected by the harness. The interface generics live in ERGM.jl; this package does `import ERGM: compute, change_stat, name` and re-exports them so user methods extend the right functions.
- **Term traits (the second half of the contract):** ERGM.jl's public term-trait protocol (`ERGM.jl/src/terms/traits.jl`) is how a term declares what it needs; ERGM's *formula validation and materialization read the traits, not term types*, so a user term is validated exactly like a built-in (this was ERGMUserterms.jl#1 — before it, only ERGM's own attribute terms were validated and user terms silently passed through). Declare what is not at its default: `ERGM.required_vertex_attributes(t) = (t.attr,)` and `ERGM.required_edge_attributes` (default `()`; declare an attribute **iff its absence is an error** — `WeightedEdges` reads `:weight` but defaults it, so it declares nothing, while `InteractionTerm` declares both of its attributes because without them it collapses to a constant zero); `ERGM.requires_directed`/`requires_undirected` (default `false`); `ERGM.is_dyad_dependent` (fallback `true`; subtyping ERGM's `NodalTerm`/`DyadicTerm` implies `false`); `Networks.supports_missing` (default `false`; `true` only if the statistic consults `is_missing_dyad` and so ignores masked dyads' face values).
- **`@ergm_term` macro:** After evaluating the body, checks that the type subtypes AbstractERGMTerm and has non-fallback compute/change_stat/name methods, warning on gaps (`_check_term_definition`).
- **Validation/testing:** `validate_term()`, `validate_traits()`, `change_stat_check()`, `consistency_check()`, `test_term()` -- verify the add-direction invariant by brute-force toggling AND check state-independence. The harness snapshots/restores the dyad's edge attributes across toggles (rem_edge! deletes them), treating them as exogenous dyadic data. `validate_traits()` (run by `validate_term` unless `traits=false`) tests the *declarations*: it perturbs each attribute (rotate values, then collapse them — rotation alone misses partition-only terms like nodematch) to check that the declared attributes are the ones the term reads (an undeclared **vertex** attribute read is a failure — that is the silent-all-zero-column trap; an undeclared **edge** attribute read is only a warning, since defaulting terms are legitimate); builds the direction twin of the network to check a `requires_directed`/`requires_undirected` declaration is enforced by `ERGMModel`; toggles *every other dyad* (not a random sample — a violation may hide in the reverse arc alone) to check an `is_dyad_dependent = false` claim; flips the face value of masked dyads to check a `supports_missing = true` claim; and finally requires `ERGMModel` construction to accept the term.
- **Package template:** `examples/MyTermPackage/` — a copyable skeleton (Project.toml + src + test) of a third-party term package: `ReciprocatedHomophily` declares all four traits, passes `validate_term`, and fits inside an `ERGMModel`. Its module is `include`d by this package's own test suite, so the template cannot rot.
- **Benchmarking:** `benchmark_term()` profiles `compute()` vs `change_stat()` performance.
- **Example terms:** `ExampleTerm`, `TemplateTerm`, `WeightedEdges`, `DyadCovTerm`, `InteractionTerm` -- shipped as reference implementations.
- **Doc helpers:** `term_signature()`, `term_documentation()` for auto-generating term docs.

## Key Dependencies

- **ERGM.jl** -- provides `AbstractERGMTerm`, core ERGM functionality
- **Networks.jl** -- network data structure (`Network`, `has_edge`, `add_edge!`, `rem_edge!`, `nv`, `ne`, etc.)
- **Graphs.jl** -- graph algorithms and iteration (`edges`, `vertices`, `src`, `dst`, `outneighbors`, `inneighbors`)
- Requires Julia >= 1.12 (Networks.jl cannot load on earlier versions)

## Conventions

- Term structs should be **immutable** (`struct`, not `mutable struct`) and subtype `AbstractUserTerm`.
- `change_stat()` should be O(degree), not O(edges), for performance.
- Term names should be descriptive and include parameters (e.g., `"template.$(t.param)"`).
- Validation uses a tolerance of `1e-10` for floating-point consistency checks.
- All exported symbols are declared at the top of the module; there are no separate include files.
