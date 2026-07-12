# Benchmarking

ERGMUserterms.jl includes performance profiling tools to ensure custom terms are efficient enough for MCMC simulation. Since `change_stat()` is called thousands of times during ERGM estimation, its performance directly affects model fitting speed.

## Why Benchmark?

In ERGM estimation via MCMC:

- Each iteration proposes toggling one edge
- `change_stat()` is called for every proposed toggle
- A typical estimation involves 10,000-1,000,000+ proposals
- Slow `change_stat()` makes estimation infeasible

The ideal `change_stat()` is much faster than `compute()`, because it only needs to calculate the *local* change rather than recomputing the full statistic.

## benchmark_term

Profile the performance of `compute()` and `change_stat()`:

<!-- skip-check -->
```julia
result = benchmark_term(term, net; n_iter=1000)
```

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
using ERGM, ERGMUserterms, Network

# Create a test network
net = network(100; directed=true)
for _ in 1:500
    i, j = rand(1:100), rand(1:100)
    i != j && add_edge!(net, i, j)
end

# Benchmark
term = ExampleTerm()
result = benchmark_term(term, net; n_iter=1000)

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

# FAST: O(1) - just checks the two vertices
function change_stat(t::NodeMatchTerm, net, i::Int, j::Int)
    attrs = get_vertex_attribute(net, t.attr)
    same = get(attrs, i, nothing) == get(attrs, j, nothing)
    return same ? 1.0 : 0.0
end
```

## Benchmarking at Different Scales: profile_term

Network size significantly affects performance. `profile_term` runs
`benchmark_term` on random networks of several sizes and returns one result
per size, which makes change statistics that do not scale (e.g. O(edges)
instead of O(degree)) easy to spot:

```julia
results = profile_term(ExampleTerm(); sizes=[10, 50, 100], density=0.1, n_iter=500)

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

net = network(50; directed=true)
for _ in 1:200
    i, j = rand(1:50), rand(1:50)
    i != j && add_edge!(net, i, j)
end

for term in terms
    result = benchmark_term(term, net; n_iter=1000)
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
result = benchmark_term(term, net; n_iter=1000)
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

## Best Practices

1. **Benchmark on realistic networks**: Use sizes similar to your actual data
2. **Aim for high speedup**: change_stat should be at least 10× faster than compute
3. **Check scaling**: Speedup should increase with network size
4. **Minimize allocations**: Avoid creating arrays or dictionaries in change_stat
5. **Use type-stable code**: Return consistent types from all code paths
6. **Profile before optimizing**: Identify the actual bottleneck first
7. **Compare with built-in terms**: Your terms should be comparable to ERGM.jl's built-in terms
8. **Re-benchmark after changes**: Performance can change with code modifications
