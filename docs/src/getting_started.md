# Getting Started

Implement a shared-neighbor statistic, verify its change calculation on several networks, and use it in a model. Keep the distinction between passing an implementation check and choosing a statistically well-posed model throughout the workflow.

!!! note "Before you begin"

    The validator checks computational consistency on the networks it exercises. It does not establish identifiability, model fit, or scientific validity. A custom term must declare its dependence and attribute requirements accurately; missing-data support is opt-in.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

Run the blocks below in order in that environment. They build on variables from earlier steps; stochastic examples use seeded random number generators where shown.

## Basic Workflow

The typical ERGMUserterms.jl workflow consists of four steps:

1. **Define the term** - Create a struct and implement the interface
2. **Validate** - Check correctness with automated tools
3. **Test** - Run comprehensive tests on various networks
4. **Use in ERGM** - Integrate with ERGM.jl for model estimation

## Step 1: Define a Custom Term

Every custom term needs a struct and three methods:

```julia
using ERGM
using ERGMUserterms
using Networks
using Graphs: src, dst                    # edge endpoints (Graphs.jl's; not re-exported by Networks)
import ERGM: name, compute, change_stat   # required to extend the term interface

# Define the struct
struct SharedNeighborTerm <: AbstractUserTerm end

# Method 1: name
name(::SharedNeighborTerm) = "shared_neighbors"

# Method 2: compute - full network statistic
function compute(::SharedNeighborTerm, net)
    total = 0.0
    for e in edges(net)
        i, j = src(e), dst(e)
        # Count common out-neighbors
        for k in outneighbors(net, i)
            k != j && has_edge(net, j, k) && (total += 1.0)
        end
    end
    return total
end

# Method 3: change_stat - the ADD-DIRECTION change statistic: the value of
# the statistic with edge (i,j) present minus with it absent, regardless of
# the dyad's current state. Neighborhoods that could contain the dyad's own
# edge are masked (k == j skips) so the value is state-independent.
function change_stat(::SharedNeighborTerm, net, i::Int, j::Int)
    delta = 0.0

    # 1. The new edge (i,j)'s own shared out-neighbors
    for k in outneighbors(net, i)
        k == j && continue
        has_edge(net, j, k) && (delta += 1.0)
    end

    # 2. Existing edges (i,b) gain shared out-neighbor j when b→j
    for b in outneighbors(net, i)
        b == j && continue
        has_edge(net, b, j) && (delta += 1.0)
    end

    # 3. Existing edges (a,i) gain shared out-neighbor j when a→j
    for a in inneighbors(net, i)
        a == j && continue
        has_edge(net, a, j) && (delta += 1.0)
    end

    return delta
end
```

### Term Struct Guidelines

| Guideline | Description |
|-----------|-------------|
| Subtype `AbstractUserTerm` | Ensures compatibility with validation tools |
| Store parameters as fields | Any configuration the term needs |
| Use type parameters for flexibility | e.g., `TemplateTerm{T}` for generic parameters |
| Keep structs immutable | Use `struct`, not `mutable struct` |

## Step 2: Validate the Term

Use `validate_term` to check your implementation:

```julia
# Create a test network
net = network(10; directed=true)
for (i, j) in [(1,2), (2,3), (1,3), (3,1), (2,1)]
    add_edge!(net, i, j)
end

# Run validation (seeded: the same dyads are checked every time, and any
# failure message ends with an rng literal that replays it)
using Random
term = SharedNeighborTerm()
valid = validate_term(term, net; verbose=true, rng=Xoshiro(1))
@assert valid    # the page's claim, checked: a rotten example fails here

# Output:
# [ Info: ✓ name() returns: shared_neighbors
# [ Info: ✓ compute() returns: 3.0 (Float64)
# [ Info: Testing change_stat() with 10 random dyads...
# [ Info: ✓ change_stat() returns valid values
# [ Info: Running consistency check...
# [ Info: ✓ change_stat() consistent with compute()
# [ Info: ✓ declared attributes match the ones the term reads
# [ Info: ✓ accepted by ERGMModel construction
```

### Validation Checks

| Check | Description |
|-------|-------------|
| `name()` | Returns a non-empty String |
| `compute()` | Returns a Real value without errors |
| `change_stat()` | Returns Real values for `n_tests` random dyads |
| Consistency | `change_stat()` matches `compute()` differences on `n_tests` random dyads |
| Traits | The term's declarations hold ([`validate_traits`](@ref)) and `ERGMModel` accepts it |

Every random dyad comes from the `rng` keyword (default `Random.default_rng()`);
pass `rng=Xoshiro(seed)` for a run you can repeat, and read the
`rng=Xoshiro(…)` literal off any failure message to replay that exact run.

## Step 3: Test Comprehensively

Use `test_term` for thorough testing on random networks. The networks it
generates are **directed and carry no vertex attributes unless you say
otherwise**: a term that declares `ERGM.required_vertex_attributes` needs
`vertex_attributes=Dict(:attr => v -> value)` (or a vector with one value per
vertex), and one that declares `ERGM.requires_undirected` needs
`directed=false` — otherwise `validate_term` reports the mismatch and the run
fails, as it should:

```julia
passed = test_term(term; n_vertices=20, density=0.1, n_tests=100, rng=Xoshiro(1))
@assert passed

# A nodematch-style term declaring `ERGM.required_vertex_attributes(t) = (t.attr,)`
# is tested on networks that carry the attribute:
@assert test_term(InteractionTerm(:a, :b); n_vertices=12, n_tests=20, rng=Xoshiro(1),
                  vertex_attributes=Dict(:a => (v -> Float64(v)), :b => (v -> Float64(13 - v))))

# Output:
# Testing term: shared_neighbors
# (replay with rng=Xoshiro(0x…, 0x…, 0x…, 0x…))
# ==================================================
# [ Info: ✓ name() returns: shared_neighbors
# [ Info: ✓ compute() returns: 6.0 (Float64)
# [ Info: Testing change_stat() with 100 random dyads...
# [ Info: ✓ change_stat() returns valid values
# [ Info: Running consistency check...
# [ Info: ✓ change_stat() consistent with compute()
# [ Info: ✓ declared attributes match the ones the term reads
# [ Info: ✓ accepted by ERGMModel construction
#
# Additional tests:
# ✓ Works on empty network: 0.0
# ✓ Works on complete network: 720.0
# ==================================================
# All tests PASSED
```

### Test Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_vertices` | Number of vertices in random test network | 20 |
| `density` | Edge density of random test network | 0.1 |
| `n_tests` | Random dyads for the `change_stat` checks (passed to `validate_term`) | 100 |
| `directed` | Build directed (`true`) or undirected (`false`) networks — `false` for a `requires_undirected` term | `true` |
| `vertex_attributes` | `attr => spec` pairs set on every generated network; `spec` is `v -> value` or a vector with one value per vertex | `Dict{Symbol,Any}()` |
| `rng` | Source of every random draw — the network and the dyads | `Random.default_rng()` |

Two things `test_term` (and every validator) will not do: validate on a
**two-mode** network — `network(n; bipartite=k)` is refused with an
`ArgumentError`, because ERGM.jl fits one-mode networks only and the
within-mode dyads cannot hold an edge — and treat a **masked dyad** as
anything but its face value for a term that does not declare
`Networks.supports_missing` (it logs one `[ Info: net has k masked dyads …]`
line so you know; ERGM.jl applies its `missing=` policy at estimation time).

## Step 4: Use in ERGM

Once validated, your term works with ERGM.jl exactly like a built-in one —
in a term list, in `ERGMModel` construction, and in `fit_ergm` (or its
statnet-style alias `ergm`, the same function):

```julia
using ERGM
using Networks
using Random

# Create a network
rng = Xoshiro(2)
net = network(50; directed=true)
for _ in 1:150
    i, j = rand(rng, 1:50), rand(rng, 1:50)
    i != j && add_edge!(net, i, j)
end

# Use your custom term next to built-ins
terms = [Edges(), Triangle(), SharedNeighborTerm()]

# Compute statistics
for term in terms
    println(name(term), ": ", compute(term, net))
end

# Fit a model containing it (MPLE by default; pass method=:mcmle for a
# Monte-Carlo MLE, and always an rng)
fit = fit_ergm(net, [Edges(), SharedNeighborTerm()]; rng=Xoshiro(3))
coef(fit)                 # two coefficients: edges, shared_neighbors
@assert length(coef(fit)) == 2
@assert ergm === fit_ergm
```

## Complete Example

```julia
using ERGM
using ERGMUserterms
using Networks
using Random
import ERGM: name, compute, change_stat   # required to extend the term interface

# === Define a parameterized term ===
struct WeightedDensity <: AbstractUserTerm
    weight::Float64
end

name(t::WeightedDensity) = "wdensity.$(t.weight)"

function compute(t::WeightedDensity, net)
    n = nv(net)
    max_edges = n * (n - 1)
    return max_edges > 0 ? t.weight * ne(net) / max_edges : 0.0
end

function change_stat(t::WeightedDensity, net, i::Int, j::Int)
    # Add-direction: one extra edge adds weight/max_edges
    n = nv(net)
    max_edges = n * (n - 1)
    max_edges == 0 && return 0.0
    return t.weight / max_edges
end

# === Validate ===
rng = Xoshiro(1)
net = network(20; directed=true)
for _ in 1:38
    i, j = rand(rng, 1:20), rand(rng, 1:20)
    i != j && add_edge!(net, i, j)
end

term = WeightedDensity(10.0)
valid = validate_term(term, net; verbose=true, rng=rng)
@assert valid

# === Run comprehensive tests ===
@assert test_term(term; n_vertices=30, density=0.15, rng=rng)

# === Benchmark performance ===
result = benchmark_term(term, net; n_iter=1000, rng=rng)
println("compute() mean: ", round(result.compute_mean * 1e6, digits=2), " μs")
println("change_stat() mean: ", round(result.change_stat_mean * 1e6, digits=2), " μs")
println("Speedup: ", round(result.speedup, digits=1), "×")

# === Use in model ===
if valid
    println("\nTerm ready for ERGM estimation!")
    println("Statistic value: ", compute(term, net))
end
```

## Working with the @ergm_term Macro

The `@ergm_term` macro checks, right after the definition, that the type
subtypes `AbstractUserTerm` and has its own `name`/`compute`/`change_stat`
methods (a bare `using` without the `import ERGM: …` line above would leave
local functions behind, and the macro warns about each missing method):

```julia
@ergm_term MyDegreeVar begin
    struct MyDegreeVar <: AbstractUserTerm end

    name(::MyDegreeVar) = "degree_variance"

    function compute(::MyDegreeVar, net)
        n = nv(net)
        n == 0 && return 0.0
        degrees = [length(outneighbors(net, v)) for v in vertices(net)]
        mean_deg = sum(degrees) / n
        return sum((d - mean_deg)^2 for d in degrees) / n
    end

    function change_stat(::MyDegreeVar, net, i::Int, j::Int)
        # Add-direction by brute force: statistic with edge (i,j) forced
        # present minus forced absent, then restore the original state.
        # (Fine for a small demo; real terms should compute the delta
        # directly in O(degree) — see the efficiency guidelines.)
        had = has_edge(net, i, j)
        had && rem_edge!(net, i, j)
        without = compute(MyDegreeVar(), net)
        add_edge!(net, i, j)
        with = compute(MyDegreeVar(), net)
        had || rem_edge!(net, i, j)
        return with - without
    end
end

@assert validate_term(MyDegreeVar(), net; verbose=false, rng=Xoshiro(4))
```

## Best Practices

1. **Start with a template**: Copy an example term and modify it
2. **Validate early**: Run `validate_term` after implementing each method
3. **Test edge cases**: Verify on empty and complete networks
4. **Optimize change_stat**: It should be O(degree), not O(edges)
5. **Use descriptive names**: Include parameters in the name string
6. **Document your terms**: Use `term_documentation()` for auto-generated docs

## Next Steps

- Learn the [Term Interface](guide/term_interface.md) in detail
- Study the [Templates and Examples](guide/templates.md) provided
- Master [Validation and Testing](guide/validation.md) tools
- Profile with [Benchmarking](guide/benchmarking.md) utilities
