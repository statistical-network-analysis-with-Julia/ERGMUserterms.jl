# Benchmarking

ERGMUserterms.jl includes timing tools (`benchmark_term`, `profile_term`) to ensure custom terms are efficient enough for MCMC simulation. Since `change_stat()` is called thousands of times during ERGM estimation, its performance directly affects model fitting speed.

## Why Benchmark?

In ERGM estimation via MCMC:

- Each iteration proposes toggling one edge
- `change_stat()` is called for every proposed toggle
- A typical estimation involves 10,000-1,000,000+ proposals
- Slow `change_stat()` makes estimation infeasible

The ideal `change_stat()` is much faster than `compute()`, because it only needs to calculate the *local* change rather than recomputing the full statistic.

## benchmark_term

Time `compute()` and `change_stat()` (`benchmark_term` is a timer, not a
profiler — for a call graph use Julia's `Profile` on your term directly):

<!-- skip-check -->
```julia
result = benchmark_term(term, net; n_iter=1000, rng=Random.default_rng())
```

`n_iter` timed calls of each function; the `change_stat` calls hit random
dyads drawn from `rng`, so the same `rng` benchmarks the same dyads (the
timings themselves are, of course, whatever the machine makes of them).

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_iter` | Timed calls of each function (must be ≥ 1) | `1000` |
| `rng` | Source of the benchmarked dyads; the same `Xoshiro(seed)` times the same dyads | `Random.default_rng()` |

`benchmark_term` throws an `ArgumentError` for `n_iter < 1` or a network with
fewer than two vertices (nothing to draw a dyad from).

### Return Value

The function returns a `NamedTuple` with:

| Field | Type | Description |
|-------|------|-------------|
| `compute_mean` | `Float64` | Mean time for `compute()` (seconds) |
| `compute_std` | `Float64` | Standard deviation of `compute()` time |
| `change_stat_mean` | `Float64` | Mean time for `change_stat()` (seconds) |
| `change_stat_std` | `Float64` | Standard deviation of `change_stat()` time |
| `speedup` | `Float64` | Ratio: `compute_mean / change_stat_mean` |

### Basic Usage

```julia
using ERGM, ERGMUserterms, Networks, Random

# Create a test network
rng = Xoshiro(1)
net = network(100; directed=true)
for _ in 1:500
    i, j = rand(rng, 1:100), rand(rng, 1:100)
    i != j && add_edge!(net, i, j)
end

# Benchmark
term = ExampleTerm()
result = benchmark_term(term, net; n_iter=1000, rng=rng)

println("compute() mean:     $(round(result.compute_mean * 1e6, digits=2)) μs")
println("change_stat() mean: $(round(result.change_stat_mean * 1e6, digits=2)) μs")
println("Speedup:            $(round(result.speedup, digits=1))×")
```

### Interpreting Results

| Speedup | Assessment | Action |
|---------|------------|--------|
| > 100× | Excellent | Term is well-optimized |
| 10-100× | Good | Acceptable for most networks |
| 2-10× | Moderate | Consider optimizing for large networks |
| ~1× | Poor | change_stat likely recomputes the full statistic |
| < 1× | Very poor | change_stat is slower than compute (bug?) |

## Performance Guidelines

### Target Complexity

| Function | Target Complexity | Description |
|----------|-------------------|-------------|
| `compute()` | O(m) or O(m + n) | Iterate over edges or vertices once |
| `change_stat()` | O(d) or O(d²) | Only examine neighbors of i and j |

Where m = number of edges, n = number of vertices, d = average degree.

### Optimizing change_stat()

The key insight: when toggling edge `(i,j)`, only statistics involving vertices `i` or `j` change. You only need to compute the *local* effect.

#### Example: Counting Triangles

<!-- skip-check -->
```julia
# SLOW: O(m) - recomputes full statistic
function change_stat(::TriangleTerm, net, i::Int, j::Int)
    # force edge absent, compute; force present, compute; restore
    # ...
    return with_edge - without_edge   # add-direction, but O(m)
end

# FAST: O(d) - only counts shared neighbors of i and j
# (add-direction: no sign flip based on the dyad's current state)
function change_stat(::TriangleTerm, net, i::Int, j::Int)
    shared = 0
    for k in outneighbors(net, i)
        k != j && has_edge(net, j, k) && (shared += 1)
    end
    return Float64(shared)
end
```

#### Example: Node Match

<!-- skip-check -->
```julia
# SLOW: O(m) - iterates all edges
function change_stat(t::NodeMatchTerm, net, i::Int, j::Int)
    before = compute(t, net)
    # ...
    return after - before
end

# SLOW in a different way: O(1) lookups, but the whole-Dict getter
# allocates an empty default Dict on every call and infers Any
# (see the pins below) — do not copy this
function change_stat(t::NodeMatchTerm, net, i::Int, j::Int)
    attrs = get_vertex_attribute(net, t.attr)
    same = get(attrs, i, nothing) == get(attrs, j, nothing)
    return same ? 1.0 : 0.0
end

# FAST: O(1), 0 B - the per-vertex getter returns `nothing` when absent
# and allocates nothing; the result is a Bool, so inference closes
function change_stat(t::NodeMatchTerm, net, i::Int, j::Int)
    a = get_vertex_attribute(net, t.attr, i)
    b = get_vertex_attribute(net, t.attr, j)
    return (a !== nothing && a == b) ? 1.0 : 0.0
end
```

## Benchmarking at Different Scales: profile_term

Network size significantly affects performance. `profile_term` runs
`benchmark_term` on random networks of several sizes and returns one result
per size, which makes change statistics that do not scale (e.g. O(edges)
instead of O(degree)) easy to spot. The networks are built as `test_term`
builds them: pass `directed=false` for an undirected-only term and
`vertex_attributes=Dict(:attr => v -> …)` for an attribute-declaring one,
so the term is timed on the branch a real fit would run:

```julia
results = profile_term(ExampleTerm(); sizes=[10, 50, 100], density=0.1, n_iter=500,
                       rng=Xoshiro(1))

for r in results
    println("n=$(r.n_vertices), m=$(r.ne): " *
            "compute=$(round(r.compute_mean*1e6, digits=1))μs, " *
            "change_stat=$(round(r.change_stat_mean*1e6, digits=1))μs, " *
            "speedup=$(round(r.speedup, digits=1))×")
end
```

Each entry is the `benchmark_term` `NamedTuple` for that size, with
`n_vertices` (network size) and `ne` (realized edge count) added.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sizes` | Network sizes (vertex counts) to profile | `[10, 20, 40]` |
| `density` | Edge density of each random test network | `0.1` |
| `n_iter` | Iterations passed to `benchmark_term` per size | `200` |
| `rng` | Source of every random draw (the networks and the benchmarked dyads) | `Random.default_rng()` |

### Expected Scaling Behavior

For a well-implemented term:

| Network Size | compute() | change_stat() | Speedup |
|-------------|-----------|---------------|---------|
| n=10 | ~1 μs | ~0.1 μs | ~10× |
| n=100 | ~50 μs | ~1 μs | ~50× |
| n=1000 | ~5 ms | ~10 μs | ~500× |

The speedup should **increase** with network size because `compute()` scales with the network while `change_stat()` stays local.

## Comparing Terms

Benchmark multiple terms on the same network to compare:

```julia
terms = [
    ExampleTerm(),
    TemplateTerm(2.0),
    WeightedEdges(),
]

rng = Xoshiro(1)
net = network(50; directed=true)
for _ in 1:200
    i, j = rand(rng, 1:50), rand(rng, 1:50)
    i != j && add_edge!(net, i, j)
end

for term in terms
    result = benchmark_term(term, net; n_iter=1000, rng=rng)
    println("$(rpad(name(term), 20)) " *
            "compute=$(round(result.compute_mean*1e6, digits=1))μs  " *
            "change_stat=$(round(result.change_stat_mean*1e6, digits=1))μs  " *
            "speedup=$(round(result.speedup, digits=1))×")
end
```

## Profiling Tips

### Warm Up the JIT

Julia compiles functions on first call. Warm up before benchmarking:

```julia
# Warm-up calls (not timed)
compute(term, net)
change_stat(term, net, 1, 2)

# Now benchmark
result = benchmark_term(term, net; n_iter=1000, rng=Xoshiro(2))
@assert result.compute_mean > 0
```

`benchmark_term` handles this internally by using many iterations, but be aware of it for manual timing.

### Memory Allocation

Minimize allocations in `change_stat()`:

<!-- skip-check -->
```julia
# BAD: allocates on every call
function change_stat(::MyTerm, net, i::Int, j::Int)
    neighbors = collect(outneighbors(net, i))  # Allocates!
    # ...
end

# GOOD: iterates without allocation
function change_stat(::MyTerm, net, i::Int, j::Int)
    for k in outneighbors(net, i)  # No allocation
        # ...
    end
end
```

### Type Stability

Ensure your functions are type-stable:

<!-- skip-check -->
```julia
# BAD: type-unstable return
function change_stat(::MyTerm, net, i::Int, j::Int)
    if some_condition
        return 1    # Returns Int
    else
        return 0.0  # Returns Float64
    end
end

# GOOD: type-stable return
function change_stat(::MyTerm, net, i::Int, j::Int)
    if some_condition
        return 1.0  # Always Float64
    else
        return 0.0  # Always Float64
    end
end
```

## The Package's Own Benchmark Suite and Allocation Pins

`benchmark/` is a standalone BenchmarkTools environment (its `Project.toml`
sources this package as `{path = ".."}` and the siblings as
`../../ERGM.jl`, `../../Networks.jl`):

```bash
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/regression_tests.jl   # allocation gates (CI runs this)
julia --project=benchmark benchmark/benchmarks.jl         # BENCHJL rows + scaling assertion
```

`benchmarks.jl` defines `SUITE` with one `@benchmarkable` per bundled term
(and the template's `ReciprocatedHomophily`) for `compute` and for
`change_stat`, on directed Erdős–Rényi networks with the *same* mean degree
at n = 500 and n = 2000; run standalone it prints one
`BENCHJL\t<name>\t<median ns>\t<allocs>\t<bytes>` line per benchmark (what
the site's `tools/run_benchmarks.jl` consumes) and exits non-zero if any
term's per-dyad `change_stat` cost grew by more than 3× between the two
sizes — an O(edges) regression shows up as ≈4×, O(1)/O(degree) stays ≈1×.

`regression_tests.jl` (mirrored by the test suite's "Allocation regressions"
testset) pins two things for every bundled term and the template's
`ReciprocatedHomophily`, on Julia 1.12.6: **0 B** per `change_stat` and per
`compute`, and `Base.return_types(change_stat, …) == [Float64]` (likewise
`compute`).

| term | `change_stat` | how |
|:--|--:|:--|
| `ExampleTerm`, `TemplateTerm`, `DyadCovTerm` | 0 B | no attribute reads |
| `WeightedEdges` | 0 B | `get_edge_attribute(net, attr, i, j)` per edge, asserted `::Float64` |
| `InteractionTerm` | 0 B | `get_vertex_attribute(net, attr, v)` per vertex, asserted `::Float64` |
| `ReciprocatedHomophily` | 0 B | `get_vertex_attribute(net, attr, v)` per vertex |

The inference pin is the important one. The attribute stores are untyped
(`Dict{Int,Any}` / `Dict{Tuple{Int,Int},Any}`), so a term that reads them
with a bare `get(get_vertex_attribute(net, attr), v, 0.0)` infers `Any` —
which is what the shipped `WeightedEdges`/`InteractionTerm` did before
0.2.0: 96 B / 208 B of boxing per call, and, worse, ERGM's statically typed
change-statistic tuple inferred `Tuple{Float64, Any}` for *every model
containing them*, so each MPLE row and MH step went through a dynamic
convert. The per-element getters return `nothing` when the value is absent,
allocate nothing (the whole-Dict getters allocate an empty default `Dict`
on every call), and the `::Float64` assertion closes inference.

A term on a real hot path should still go one step further and snapshot its
attributes once, typed, at construction —
`vertex_attribute_vector(net, attr, Float64)` /
`get_edge_attribute(net, attr, Float64)` — and index a `Vector{Float64}` in
`change_stat`: a vector index instead of two hash lookups (see
[Attribute Validation and Snapshotting](@ref)). Any allocation, or an `Any`
return type, is a regression.

## Best Practices

1. **Benchmark on realistic networks**: Use sizes similar to your actual data
2. **Aim for high speedup**: change_stat should be at least 10× faster than compute
3. **Check scaling**: Speedup should increase with network size
4. **Minimize allocations**: Avoid creating arrays or dictionaries in change_stat
5. **Use type-stable code**: Return consistent types from all code paths
6. **Profile before optimizing**: Identify the actual bottleneck first
7. **Compare with built-in terms**: Your terms should be comparable to ERGM.jl's built-in terms
8. **Re-benchmark after changes**: Performance can change with code modifications
