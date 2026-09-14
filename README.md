# ERGMUserterms.jl


[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl)
[![Build Status](https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.svg" alt="ERGMUserterms.jl icon" width="160">
</p>

Custom ERGM Term Development for Julia.

## Overview

ERGMUserterms.jl provides templates, utilities, and validation tools for developing custom ERGM terms. It includes example terms, a comprehensive testing framework, and documentation helpers.

This package is a Julia port of the R `ergm.userterms` package from the StatNet collection.

### Not implemented (vs R ergm.userterms)

What is ported is the *role* of `ergm.userterms` — a validated, copyable
starting point for third-party terms — not its contents. R's package is a C
template: its `changestats.users.c` skeleton and its worked example term
`mindegree` are **not** ported, and nothing shipped here has them as a
counterpart (`ergm.userterms` is archived from CRAN and is not what the
bundled terms are pinned against). A Julia term is plain Julia — three
methods plus trait declarations — so there is no C skeleton to fill in; the
Julia counterpart of the C template is
[`examples/MyTermPackage/`](examples/MyTermPackage), and the bundled example
terms are pinned against the `ergm` terms they are identical to (see
[Validation against R](#validation-against-r)).

## Installation

Requires Julia 1.12+. ERGMUserterms.jl depends on the unregistered
[Networks.jl](https://github.com/statistical-network-analysis-with-Julia/Networks.jl) and [ERGM.jl](https://github.com/statistical-network-analysis-with-Julia/ERGM.jl) packages, which must be added first (in this order):

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/Networks.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGM.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl")
```

For development, you can instead clone all ecosystem repositories side by
side (the monorepo layout) and start Julia with the root workspace project
(`julia --project=.` in the clone root): the `[sources]` path dependencies
then wire the packages together with no ordered installs needed.

## Features

- **Term macro**: `@ergm_term` for defining custom terms
- **Validation**: Automatic validation of term implementations
- **Testing**: Consistency checks between `compute()` and `change_stat()`
- **Benchmarking**: Timing of `compute()` vs `change_stat()` for terms
- **Templates**: Example terms to copy and modify

## Quick Start

```julia
using ERGM
using ERGMUserterms
using Networks

# Extend the shared interface generics (required so ERGM.jl sees your
# methods; ERGMUserterms re-exports the same three functions)
import ERGM: name, compute, change_stat

# Define a custom term
struct MyTerm <: AbstractUserTerm
    param::Float64
end

name(t::MyTerm) = "myterm.$(t.param)"

function compute(t::MyTerm, net)
    # Count edges weighted by parameter
    return Float64(ne(net)) * t.param
end

function change_stat(t::MyTerm, net, i::Int, j::Int)
    # Add-direction change: statistic with edge (i,j) present minus with it
    # absent. Must NOT depend on whether the edge currently exists.
    return t.param
end

# Validate the term (seeded: the harness draws every random dyad from `rng`,
# and a failure message ends with the rng literal that replays it)
using Random
net = network(20; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
term = MyTerm(2.0)
valid = validate_term(term, net; rng=Xoshiro(1))
@assert valid
```

`import ERGM: name, compute, change_stat` and
`import ERGMUserterms: name, compute, change_stat` name the same three
functions (`ERGMUserterms.compute === ERGM.compute === Networks.compute` is
pinned by the test suite), so either spelling works — but one of them is
required: with a bare `using`, `name(::MyTerm) = …` defines a *local*
function that ERGM.jl never calls, and `validate_term` reports
`compute() failed`. Edge endpoints `src`/`dst` are `Graphs.src`/`Graphs.dst`
(not re-exported by Networks.jl): a term iterating `edges(net)` needs
`using Graphs: src, dst`.

## Term Interface

Every ERGM term must implement:

<!-- skip-check -->
```julia
# Required
name(term) -> String           # Term name for output
compute(term, net) -> Float64  # Network statistic
change_stat(term, net, i, j) -> Float64  # Add-direction change statistic
```

The key relationship (the **add-direction** convention):
```
change_stat(term, net, i, j) == compute(term, net⁺ij) - compute(term, net⁻ij)
```
where `net⁺ij`/`net⁻ij` are `net` with edge (i,j) forced present/absent and
all other dyads unchanged. The value must be the same whether or not the
edge currently exists — the toggle-direction idiom
`has_edge(net, i, j) ? -Δ : Δ` is wrong and is rejected by the validation
harness (ERGM.jl's MH sampler negates the add-direction value itself for
removal proposals). That idiom is R ergm's own C convention —
`C_CHANGESTAT_FN(…, Rboolean edgestate)`, with every statnet changestat
written `CHANGE_STAT[0] += edgestate ? -1 : 1` — so when porting an
`ergm.userterms` changestat, drop the `edgestate` sign and return the
add-direction value.

A term also **declares its traits** through ERGM.jl's public term-trait
protocol. ERGM's formula validation reads the declarations, not the term's
type, so a custom term is validated at model construction exactly like a
built-in one:

<!-- skip-check -->
```julia
ERGM.required_vertex_attributes(t::MyTerm) = (t.attr,)  # default ()
ERGM.required_edge_attributes(t::MyTerm)   = ()         # default ()
ERGM.requires_directed(::MyTerm)           = true       # default false
ERGM.requires_undirected(::MyTerm)         = false      # default false
ERGM.is_dyad_dependent(::MyTerm)           = false      # default true (conservative)
Networks.supports_missing(::MyTerm)        = true       # default false
```

Declare an attribute iff its absence is an *error*: a term reading an
undeclared attribute silently becomes an all-zero design column on a network
that lacks it, whereas a declared one raises an `ArgumentError` naming it.
Declare `is_dyad_dependent = false` only for covariate-only terms — the
fallback `true` triggers ERGM.jl's pseudo-likelihood caveat and a conservative
MCMLE bridge reference. Declare `supports_missing = true` only if the statistic
consults `is_missing_dyad` and so ignores masked dyads' face values.

Two of the traits deserve a word. A **declared vertex attribute must have a
value on every vertex**: `ERGMModel` refuses a partial one with `term '…'
needs vertex attribute :a on every vertex … statnet refuses NA` (the check is
`ERGM._validate_formula`, public in ERGM.jl), so a `get(attrs, v, default)`
fallback in your code is reached only by raw `compute` calls. And
`Networks.supports_missing` **is** `ERGM.supports_missing`
(`ERGM.supports_missing === Networks.supports_missing`): either spelling
extends the same generic, and the idiom is the one ERGM.jl's own docstring
shows — consult the mask, then declare:

```julia
using Graphs: src, dst
using Random
struct ObservedEdges <: AbstractUserTerm end
name(::ObservedEdges) = "observed_edges"
compute(::ObservedEdges, net) =
    Float64(count(!is_missing_dyad(net, src(e), dst(e)) for e in edges(net)))
change_stat(::ObservedEdges, net, i::Int, j::Int) = is_missing_dyad(net, i, j) ? 0.0 : 1.0
Networks.supports_missing(::ObservedEdges) = true   # ≡ ERGM.supports_missing(::ObservedEdges) = true

masked = network(5; directed=true)
add_edge!(masked, 1, 2); add_edge!(masked, 2, 3); set_missing_dyad!(masked, 1, 2)
compute(ObservedEdges(), masked)                    # 1.0 — the masked tie does not count
@assert validate_traits(ObservedEdges(), masked; verbose=false, rng=Xoshiro(1))
```

`validate_term` exercises all of it, and
[`examples/MyTermPackage/`](examples/MyTermPackage) is a copyable package
template for a third-party term declaring the lot.

## Validation

<!-- skip-check -->
```julia
# Full validation: name(), compute(), change_stat() on n_tests random dyads,
# consistency, and the trait declarations (validate_traits)
validate_term(term, net; verbose=true, n_tests=10, rng=Xoshiro(1))

# Just consistency check
change_stat_check(term, net; n_tests=10, rng=Xoshiro(1))

# Exhaustive consistency (slow for large networks)
consistency_check(term, net; exhaustive=true)
```

Every validator, tester and benchmark takes `rng::AbstractRNG` (default
`Random.default_rng()`) and draws nothing from anywhere else: two calls with
equal rngs check the same dyads and return the same verdict. A failure
message ends with the seed line — the rng state on entry — and passing it
back replays that exact run, failing dyad included:

<!-- skip-check -->
```julia
change_stat_check(broken_term, net; verbose=true)
# ┌ Warning: State-dependent change_stat at dyad (2,5): -1.0 with edge state
# │ as-is vs 1.0 after toggling. change_stat must return the add-direction
# │ change regardless of whether the edge exists; reproduce with
# │ rng=Xoshiro(0x62bea62705299f9d, 0xb2196eb285598f6a, 0x813136b07248b437, 0x783171c16f41f41d)
# │ passed to change_stat_check
# false
```

| keyword | on | default | meaning |
|:--|:--|:--|:--|
| `rng` | every validator, `test_term`, `benchmark_term`, `profile_term` | `Random.default_rng()` | source of every random draw |
| `n_tests` | `validate_term`, `change_stat_check`, `test_term` | 10 / 10 / 100 | random dyads for the `change_stat` checks (at least 1; `test_term` now honours it; capped at the number of dyads) |
| `verbose` | validators | `true` (`consistency_check`: `false`) | log each check and failure; the returned flag is the same either way |
| `directed`, `vertex_attributes` | `test_term`, `profile_term` | `true`, `Dict{Symbol,Any}()` | the kind of network they generate: directedness, and `attr => v -> value` / `attr => vector` attributes set on every generated network |

Three limits. **No verdict without evidence**: `n_tests` must be at least 1,
and a network with fewer than 2 vertices (no dyad to check) is refused with
an `ArgumentError` by every validator and by `benchmark_term` — as are
`test_term(…; n_vertices=1)` and a `density` outside `[0, 1]`; a run that
checks nothing returns no PASS. **Networks ERGM.jl cannot fit are refused up
front**, with ERGM.jl's own message, never validated dyad by dyad and then
rejected in words that blame the term: a two-mode network
(`network(n; bipartite=k)` — the brute-force reference would read 0 on the
within-mode dyads, which cannot hold an edge, and report a correct term
inconsistent) and a network containing a **self-loop** (ERGM.jl's
statistics would count it while its estimators and samplers never touch the
diagonal; `rem_edge!(net, v, v)` is the fix the message names). **Masked
dyads are face value** for
a term that does not declare `Networks.supports_missing`: the checks run on
what the term computes, and with `verbose=true` one `[ Info: net has k
masked dyads; the term declares supports_missing = false, so they are
validated at their face value …]` line says so (ERGM.jl applies its
`missing=` policy at estimation time).

### Validation against R

The bundled terms are what you copy, so they are pinned against statnet
rather than against hand-typed literals: `test/fixtures/userterms_examples.toml`
(generated by the checked-in `test/fixtures/r/userterms_examples.R` with
ergm 4.12.0 on R 4.6.1, `[provenance]` block included) records `summary()`
and the full `ergmMPLE(output="array")` change-statistic array of the
statnet terms each bundled term is mathematically identical to, on an
8-vertex directed network and on its undirected projection:

| ERGMUserterms.jl term          | statnet term                        |
|:-------------------------------|:------------------------------------|
| `ExampleTerm()`                | `nodecov("id")` with `id` = vertex index |
| `TemplateTerm(1.0)`            | `edges`                             |
| `WeightedEdges()`              | `edgecov(W)` (2.5 on the weighted arcs, 1.0 = the default elsewhere) |
| `DyadCovTerm(cov)`             | `edgecov(cov)`, asymmetric `cov[i,j] = 2i + j` |
| `InteractionTerm(:a, :b)`      | `edgecov(M)`, `M[i,j] = a_i b_j + a_j b_i` |
| `ReciprocatedHomophily(:group)` (the `examples/MyTermPackage` template) | `mutual(same="group", diff=FALSE)` |

The test suite compares every `compute` and every `change_stat` — on every
ordered dyad, in both dyad orders on the undirected network — **and the
harness's own brute-force reference** (`_brute_change_stat`, the number
`validate_term` holds your term to) against those R values at `1e-9`, so the
harness is validated against ergm's C change statistics rather than trusted.
R's coefficient *labels* (`nodecov.id`, `edgecov.W`, …) are recorded but not
claimed: user terms keep their own names.

## Testing

```julia
# Comprehensive test suite: a random network (from rng), validate_term on it
# with n_tests random dyads, then empty and complete networks
passed = test_term(term; n_vertices=20, density=0.1, n_tests=100, rng=Xoshiro(1))
@assert passed

# The generated networks are directed and attribute-free unless you say
# otherwise: a term declaring vertex attributes gets them through
# `vertex_attributes` (a function of the vertex index, or a vector), a
# `requires_undirected` term gets `directed=false`
@assert test_term(InteractionTerm(:a, :b); n_vertices=12, n_tests=20, rng=Xoshiro(1),
                  vertex_attributes=Dict(:a => (v -> Float64(v)), :b => (v -> Float64(13 - v))))
@assert test_term(TemplateTerm(2.0); n_vertices=12, n_tests=20, directed=false, rng=Xoshiro(1))
```

## Benchmarking

```julia
result = benchmark_term(term, net; n_iter=1000, rng=Xoshiro(1))
# Returns (seconds):
#   compute_mean, compute_std
#   change_stat_mean, change_stat_std
#   speedup (compute/change_stat ratio)
@assert result.compute_mean > 0
```

`benchmark_term` times; it does not profile. `profile_term(term; sizes=…,
rng=…)` runs it on random networks of several sizes so an O(edges)
`change_stat` shows up as a cost that grows with `n`. It builds those
networks exactly as `test_term` does — `directed=false` for an
undirected-only term, `vertex_attributes=Dict(:attr => v -> …)` for an
attribute-declaring one, so the term is timed on its real branch rather
than on an attribute-absent fallback.

The package's own `benchmark/` directory is a BenchmarkTools suite over the
bundled terms (`compute` and `change_stat` at n = 500 and 2000 with constant
mean degree, asserting the per-dyad cost does not grow with n) plus
`regression_tests.jl`, the allocation pins CI runs:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/regression_tests.jl   # allocation gates
julia --project=benchmark benchmark/benchmarks.jl         # BENCHJL rows + scaling
```

Every bundled term and the template's `ReciprocatedHomophily` is pinned at
**0 B** per `change_stat` and per `compute`, and pinned to **infer
`Float64`** — the attribute stores are untyped `Dict{…,Any}`s, and a term
that reads them with a bare `get` infers `Any`, which boxes on every call
and poisons ERGM's statically typed change-statistic tuple for the whole
model. The bundled terms read per element (`get_vertex_attribute(net, attr,
v)` / `get_edge_attribute(net, attr, i, j)`, which return `nothing` when
absent and allocate nothing) and assert `::Float64`. A real term on a hot
path should go one step further and snapshot its attributes once, typed, at
construction (`vertex_attribute_vector(net, attr, Float64)` /
`get_edge_attribute(net, attr, Float64)`) — see the
[term interface guide](https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl/dev/guide/term_interface/).

## Example Terms

### ExampleTerm
<!-- skip-check -->
```julia
# Edges weighted by vertex ID sum
struct ExampleTerm <: AbstractUserTerm end

compute(::ExampleTerm, net) = sum(src(e) + dst(e) for e in edges(net))
change_stat(::ExampleTerm, net, i, j) = Float64(i + j)  # add-direction, state-independent
```

### TemplateTerm
<!-- skip-check -->
```julia
# Parameterized template
struct TemplateTerm{T} <: AbstractUserTerm
    param::T
    attr::Symbol
end

# Copy and modify for your own terms
```

### WeightedEdges
<!-- skip-check -->
```julia
# Sum of edge weights, `default` for an edge without one — read PER EDGE
# with the allocation-free getter (the whole-Dict getter infers Any)
struct WeightedEdges <: AbstractUserTerm
    attr::Symbol
    default::Float64
end

function _edge_weight(t::WeightedEdges, net, i::Int, j::Int)
    w = get_edge_attribute(net, t.attr, i, j)      # nothing when absent
    return w === nothing ? t.default : Float64(w)::Float64
end
compute(t::WeightedEdges, net) =
    sum(_edge_weight(t, net, Int(src(e)), Int(dst(e))) for e in edges(net); init=0.0)
change_stat(t::WeightedEdges, net, i::Int, j::Int) = _edge_weight(t, net, i, j)
```

### DyadCovTerm
<!-- skip-check -->
```julia
# Dyadic covariate
struct DyadCovTerm <: AbstractUserTerm
    covariate::Matrix{Float64}
end

compute(t::DyadCovTerm, net) = sum(t.covariate[src(e), dst(e)] for e in edges(net))
```

### InteractionTerm
<!-- skip-check -->
```julia
# Interaction between two node attributes
struct InteractionTerm <: AbstractUserTerm
    attr1::Symbol
    attr2::Symbol
end
```

## Documentation Helpers

```julia
# Generate term signature
sig = term_signature(term)
# "MyTerm(param::Float64)"

# Generate documentation
doc = term_documentation(term)
```

## Best Practices

1. **Efficiency**: `change_stat()` should be O(degree) not O(edges)
2. **Consistency**: Always verify with `change_stat_check()`
3. **Edge cases**: Test on empty and complete networks
4. **Naming**: Use descriptive names with parameters
5. **Reproducibility**: Pass `rng=Xoshiro(seed)` in your test suite, so a red run is the same run every time

## Common Patterns

### Counting Subgraphs
<!-- skip-check -->
```julia
function compute(::TriangleTerm, net)
    count = 0
    for i in vertices(net)
        for j in neighbors(net, i)
            for k in neighbors(net, j)
                k > i && has_edge(net, i, k) && (count += 1)
            end
        end
    end
    return count / 3  # Each triangle counted 3 times
end
```

### Using Attributes
<!-- skip-check -->
```julia
# Read per vertex: `get_vertex_attribute(net, attr, v)` returns `nothing`
# when the vertex has no value and allocates nothing; the whole-Dict
# `get_vertex_attribute(net, attr)` would allocate and infer Any (see
# term_interface.md, "Attribute Validation and Snapshotting")
function compute(t::NodeMatchTerm, net)
    count = 0.0
    for e in edges(net)
        a = get_vertex_attribute(net, t.attr, Int(src(e)))
        b = get_vertex_attribute(net, t.attr, Int(dst(e)))
        a !== nothing && a == b && (count += 1.0)
    end
    return count
end
ERGM.required_vertex_attributes(t::NodeMatchTerm) = (t.attr,)   # absent => error, not zeros
```

## Documentation

For more detailed documentation, see:

- [Stable Documentation](https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl/stable/)
- [Development Documentation](https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl/dev/)

## References

1. Hunter, D.R., Goodreau, S.M., Handcock, M.S. (2013). ergm.userterms: A template package for extending statnet. *Journal of Statistical Software*, 52(2), 1-25.

2. Hunter, D.R., Handcock, M.S., Butts, C.T., Goodreau, S.M., Morris, M. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. *Journal of Statistical Software*, 24(3), 1-29.

## Citation

If you use ERGMUserterms.jl in your work, please cite it using the entry in
[`CITATION.bib`](CITATION.bib):

```biblatex
@misc{SNWJERGMUsertermsJL,
  author = {{Statistical Network Analysis with Julia}},
  title = {ERGMUserterms.jl: Custom ERGM Term Development for Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.
