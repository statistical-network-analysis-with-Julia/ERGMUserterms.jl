# ERGMUserterms.jl


[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl)
[![Build Status](https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl/actions/workflows/CI.yml?query=branch%3Amain)
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

Limits of what the validation harness establishes:

- **A pass is a statement about the network you validated on.** By default
  every dyad of that network is checked (and the path to it from the empty
  network), but a change statistic that is wrong only in a configuration the
  network does not contain is not seen. Validate on networks that contain
  the configurations your term counts, and use `test_term` for random ones.
- **Above 5000 dyads the per-dyad checks sample** (stratified by dyad state,
  5000 dyads), and the dyad-independence check toggles every pair of dyads
  only up to about 500 dyads (22 vertices directed, 32 undirected); beyond
  that it toggles the reverse arc and the dyads sharing an endpoint first and
  a random remainder. The log line states how many dyads were checked.
- **Edge attributes are not preserved by ERGM.jl's samplers.** A term that
  reads one must hold it as a matrix (as `WeightedEdges` and `DyadCovTerm`
  do). `validate_traits` fails a term that reads it live, and ERGM.jl's
  samplers (`mcmle`, `simulate_ergm`, `gof`, the MPLE bootstrap) refuse one
  with an `ArgumentError`; the MPLE, which never toggles a tie, accepts it.
  ERGM.jl does trust a term's `is_dyad_dependent` declaration beyond a
  probe at model construction, which is why the validator tests that claim.
- The harness checks computation, not statistics: it does not establish
  identifiability, non-degeneracy or model fit.

## Installation

Requires Julia 1.12 or newer. The packages are not yet registered.

**Recommended: the ecosystem workspace.** It clones every package side by
side, develops them together in one environment, and adds the packages the
examples also use (CSV, DataFrames, Distributions, Graphs, StatsAPI,
StatsBase):

```bash
mkdir network-analysis && cd network-analysis
git clone https://github.com/statistical-network-analysis-with-Julia/statistical-network-analysis-with-Julia.github.io
julia statistical-network-analysis-with-Julia.github.io/tools/prepare_workspace.jl "$PWD" --clone
julia --project=.snippet-env
```

**Only this package, in your own environment.** Add its dependencies first,
in this order, then the extra the examples below use (`Graphs`, for
`src`/`dst`):

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/NetworkCore.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGM.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl")
Pkg.add("Graphs")
```

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
using NetworkCore

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

# Validate the term: every dyad of the network is checked against a
# brute-force recomputation (whatever is random is drawn from `rng`, and a
# failure message ends with the rng literal that replays it)
using Random
net = network(20; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
term = MyTerm(2.0)
valid = validate_term(term, net; rng=Xoshiro(1))
@assert valid

# ... and use it like any built-in term. It is an edge count, so it can
# say it is dyad-independent; the default method=:auto then fits the MPLE,
# which is the exact MLE (an undeclared term counts as dyad-dependent and
# gets the MCMLE, as R's ergm() would)
ERGM.is_dyad_dependent(::MyTerm) = false
fit = fit_ergm(load_dataset(:florentine_marriage), [MyTerm(1.0)])
@assert fit.converged && fit.method === :mple
```

`import ERGM: name, compute, change_stat` and
`import ERGMUserterms: name, compute, change_stat` name the same three
functions (`ERGMUserterms.compute === ERGM.compute === NetworkCore.compute` is
pinned by the test suite), so either spelling works — but one of them is
required: with a bare `using`, `name(::MyTerm) = …` defines a *local*
function that ERGM.jl never calls, and `validate_term` reports
`compute() failed`. Edge endpoints `src`/`dst` are `Graphs.src`/`Graphs.dst`
(not re-exported by NetworkCore.jl): a term iterating `edges(net)` needs
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
NetworkCore.supports_missing(::MyTerm)        = true       # default false
```

Declare an attribute iff its absence is an *error*: a term reading an
undeclared attribute silently becomes an all-zero design column on a network
that lacks it, whereas a declared one raises an `ArgumentError` naming it.
Declare `is_dyad_dependent = false` only for covariate-only terms — the
fallback `true` triggers ERGM.jl's pseudo-likelihood caveat and a conservative
MCMLE bridge reference. ERGM.jl trusts the declaration (`is_exact`, the exact
log-normaliser), so `validate_traits` tests it hard: a term declared
dyad-independent whose change statistic moves with any other dyad fails. Declare `supports_missing = true` only if the statistic
consults `is_missing_dyad` and so ignores masked dyads' face values.

Two of the traits deserve a word. A **declared vertex attribute must have a
value on every vertex**: `ERGMModel` refuses a partial one with `term '…'
needs vertex attribute :a on every vertex … statnet refuses NA` (the check is
`ERGM.Extension.validate_formula`, ERGM.jl's extension API), so a `get(attrs, v, default)`
fallback in your code is reached only by raw `compute` calls. And
`NetworkCore.supports_missing` **is** `ERGM.supports_missing`
(`ERGM.supports_missing === NetworkCore.supports_missing`): either spelling
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
NetworkCore.supports_missing(::ObservedEdges) = true   # ≡ ERGM.supports_missing(::ObservedEdges) = true

masked = network(5; directed=true)
add_edge!(masked, 1, 2); add_edge!(masked, 2, 3); set_missing_dyad!(masked, 1, 2)
compute(ObservedEdges(), masked)                    # 1.0 — the masked tie does not count
@assert validate_traits(ObservedEdges(), masked; verbose=false, rng=Xoshiro(1))
```

`validate_term` exercises all of it, and
[`examples/MyTermPackage/`](examples/MyTermPackage) is a copyable package
template for a third-party term declaring the lot.

## Validation

```julia
# Full validation: name(), compute(), change_stat() on EVERY dyad against a
# brute-force recomputation, g(∅) and the path from the empty network, and
# the trait declarations (validate_traits)
@assert validate_term(term, net; verbose=false)

# Just the per-dyad consistency check
@assert change_stat_check(term, net)

# A cheaper, weaker run on a big network: a stratified sample of dyads
@assert validate_term(term, net; verbose=false, n_tests=50, rng=Xoshiro(1))

# The same comparison, stopping at the first inconsistent dyad
@assert consistency_check(term, net)
```

What `validate_term` holds a term to:

- **Every dyad, by default.** A change statistic is usually wrong on a
  *kind* of dyad (the reciprocated ones, the ones closing a triangle), a
  small minority in a sparse network. `n_tests=nothing` checks them all (a
  stratified 5000 on a network with more dyads than that); an integer
  `n_tests` is a stratified sample — by the dyad's state, the reverse arc's
  state and whether the endpoints share a partner — and the log says how
  many dyads were checked.
- **`g(∅)` and the path.** `compute` on the empty network plus the change
  statistics summed while the edges are added one at a time must equal
  `compute` on the network. `g(∅)` need not be 0 (an isolates count is `n`
  there; ERGM.jl's likelihood includes `θ'g(∅)`), but a
  `ne(net) == 0 && return 0.0` shortcut for such a statistic fails.
- **The declarations**, including that an `is_dyad_dependent = false` claim
  holds on the empty, the observed and the complete network and under
  toggling of the other dyads, and that **edge attributes survive the
  sampler**: `rem_edge!` deletes a dyad's edge attributes and ERGM.jl's
  samplers toggle with it, so a term reading an edge attribute live would be
  wrong under MCMLE, `simulate_ergm` and `gof` (ERGM.jl's samplers refuse
  such a term with an `ArgumentError`, so the harness says so before a fit
  does). Hold the attribute
  as a matrix — in the constructor (`DyadCovTerm`) or through
  `ERGM.Extension.materialize(term, net)` (`WeightedEdges`).

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
| `n_tests` | `validate_term`, `change_stat_check`, `test_term` | `nothing` | `nothing`: every dyad (a stratified 5000 beyond 5000 dyads); an integer (at least 1): a stratified sample of that many distinct dyads |
| `exhaustive` | `consistency_check` | `true` | `false`: up to 100 random dyads instead of all |
| `verbose` | validators | `true` (`consistency_check`: `false`) | log each check and failure; the returned flag is the same either way |
| `directed`, `vertex_attributes` | `test_term`, `profile_term` | `true`, `Dict{Symbol,Any}()` | the kind of network they generate: directedness, and `attr => v -> value` / `attr => vector` attributes set on every generated network |

Three limits. **No verdict without evidence**: an integer `n_tests` must be at least 1,
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
a term that does not declare `NetworkCore.supports_missing`: the checks run on
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
claimed: user terms keep their own names. The fixture also records the
**fitted** `edgecov(W)` model on both networks; `WeightedEdges()` reproduces
R's coefficient, standard error and log-likelihood by MPLE and its
coefficient by MCMLE, and its sampler draws are, draw for draw, those of
ERGM's `EdgeCov(W)`.

## Testing

```julia
# Comprehensive test suite: a random network (from rng), validate_term on
# every dyad of it, then empty and complete networks
passed = test_term(term; n_vertices=20, density=0.1, rng=Xoshiro(1))
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

The bundled `ExampleTerm` (statnet's `nodecov("id")`), written out under
another name so the block runs as it stands — defined, validated and fitted:

```julia
using Graphs: src, dst

# Edges weighted by vertex ID sum
struct IdSum <: AbstractUserTerm end

name(::IdSum) = "idsum"
compute(::IdSum, net) = sum(Float64(src(e) + dst(e)) for e in edges(net); init=0.0)
change_stat(::IdSum, net, i, j) = Float64(i + j)  # add-direction, state-independent
ERGM.is_dyad_dependent(::IdSum) = false           # covariate-only

flo = load_dataset(:florentine_marriage)
@assert validate_term(IdSum(), flo; verbose=false)
fit = fit_ergm(flo, [Edges(), IdSum()])
@assert coef(fit) ≈ coef(fit_ergm(flo, [Edges(), ExampleTerm()]))
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

Sum of edge weights, `default` for a dyad without one: statnet's
`edgecov(W)`. Inside a model the weights are a matrix snapshotted when the
model is built, because the samplers' toggles delete live edge attributes:

```julia
wnet = network(5; directed=true)
add_edge!(wnet, 1, 2); add_edge!(wnet, 2, 3)
set_edge_attribute!(wnet, :weight, 1, 2, 3.0)
w = WeightedEdges()                                # attr = :weight, default = 1.0
@assert compute(w, wnet) == 4.0                    # 3.0 + the default 1.0
held = ERGMModel(ERGMFormula([w]), wnet).formula.terms[1]   # what a fit evaluates
rem_edge!(wnet, 1, 2); add_edge!(wnet, 1, 2)       # a sampler toggle: the stored weight is gone
@assert change_stat(held, wnet, 1, 2) == 3.0       # the model still sees edgecov(W)
```

The [templates guide](https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl/dev/guide/templates/)
shows the full source, including the `ERGM.Extension.materialize` method that takes
the snapshot.

### DyadCovTerm
<!-- skip-check -->
```julia
# Dyadic covariate: an n×n matrix in vertex order. A matrix of another size
# is an ArgumentError from compute, change_stat and ERGMModel construction.
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
3. **Edge cases**: Test on empty and complete networks; never special-case the empty network to 0 unless the statistic is 0 there
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

- [Development Documentation](https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl/dev/)

## References

1. Hunter, D.R., Goodreau, S.M., Handcock, M.S. (2013). ergm.userterms: A template package for extending statnet. *Journal of Statistical Software*, 52(2), 1-25.

2. Hunter, D.R., Handcock, M.S., Butts, C.T., Goodreau, S.M., Morris, M. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. *Journal of Statistical Software*, 24(3), 1-29.

3. Krivitsky, P.N., Hunter, D.R., Morris, M., Klumb, C. (2023). ergm 4: New features for analyzing exponential-family random graph models. *Journal of Statistical Software*, 105(6), 1-44.

## Citation

If you use ERGMUserterms.jl in your work, please cite it using the entry in
[`CITATION.bib`](CITATION.bib). Please also cite the R packages it is a
counterpart of and their methods papers — `ergm.userterms` (Hunter,
Goodreau & Handcock 2013) and `ergm` (Hunter et al. 2008; Krivitsky et al.
2023); the ecosystem's
[How to cite](https://statistical-network-analysis-with-julia.github.io/citing/)
page lists the references with DOIs.

```biblatex
@misc{SNWJERGMUsertermsJL,
  author = {Santoni, Simone},
  title = {ERGMUserterms.jl: Custom ERGM Term Development for Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.
