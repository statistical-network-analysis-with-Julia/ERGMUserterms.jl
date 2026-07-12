# Getting Started

This tutorial walks through creating, validating, and using custom ERGM terms with ERGMUserterms.jl.

## Installation

Install ERGMUserterms.jl from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/Network.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGM.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl")
```

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
using Network
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

# Run validation
term = SharedNeighborTerm()
valid = validate_term(term, net; verbose=true)

# Output:
# [ Info: ✓ name() returns: shared_neighbors
# [ Info: ✓ compute() returns: 4.0 (Float64)
# [ Info: Testing change_stat() with 10 random edges...
# [ Info: ✓ change_stat() returns valid values
# [ Info: Running consistency check...
# [ Info: ✓ change_stat() consistent with compute()
```

### Validation Checks

| Check | Description |
|-------|-------------|
| `name()` | Returns a non-empty String |
| `compute()` | Returns a Real value without errors |
| `change_stat()` | Returns Real values for random edges |
| Consistency | `change_stat()` matches `compute()` differences |

## Step 3: Test Comprehensively

Use `test_term` for thorough testing on random networks:

```julia
test_term(term; n_vertices=20, density=0.1, n_tests=100)

# Output:
# Testing term: shared_neighbors
# ==================================================
# [ Info: ✓ name() returns: shared_neighbors
# [ Info: ✓ compute() returns: 12.0 (Float64)
# [ Info: Testing change_stat() with 10 random edges...
# [ Info: ✓ change_stat() returns valid values
# [ Info: Running consistency check...
# [ Info: ✓ change_stat() consistent with compute()
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
| `n_tests` | Number of random edges for consistency | 100 |

## Step 4: Use in ERGM

Once validated, your term works with ERGM.jl:

```julia
using ERGM
using Network

# Create a network
net = network(50; directed=true)
# ... add edges ...

# Use your custom term in a model
terms = [Edges(), Triangle(), SharedNeighborTerm()]

# Compute statistics
for term in terms
    println(name(term), ": ", compute(term, net))
end
```

## Complete Example

```julia
using ERGM
using ERGMUserterms
using Network
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
net = network(20; directed=true)
for _ in 1:38
    i, j = rand(1:20), rand(1:20)
    i != j && add_edge!(net, i, j)
end

term = WeightedDensity(10.0)
valid = validate_term(term, net; verbose=true)

# === Run comprehensive tests ===
test_term(term; n_vertices=30, density=0.15)

# === Benchmark performance ===
result = benchmark_term(term, net; n_iter=1000)
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

The `@ergm_term` macro provides a convenient way to define terms with automatic validation:

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
