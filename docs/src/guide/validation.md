# Validation and Testing

ERGMUserterms.jl provides a comprehensive validation framework to ensure custom terms are correctly implemented before use in ERGM estimation. Incorrect terms lead to biased estimates, so validation is critical.

## Overview

The validation process has three levels:

1. **Quick validation**: `validate_term()` checks all methods work and are consistent
2. **Consistency checking**: `change_stat_check()` and `consistency_check()` verify the fundamental relationship
3. **Comprehensive testing**: `test_term()` runs all checks on multiple network types

## validate_term

The primary validation function that checks all aspects of a term implementation:

```julia
valid = validate_term(term, net; verbose=true)
```

### What It Checks

| Check | Description | Failure Means |
|-------|-------------|---------------|
| `name()` exists | Returns a non-empty string | Method not implemented |
| `name()` non-empty | String length > 0 | Empty name will cause display issues |
| `compute()` works | Returns without error | Bug in compute implementation |
| `compute()` returns Real | Type check on return value | Wrong return type |
| `change_stat()` works | Returns without error for random edges | Bug in change_stat implementation |
| `change_stat()` returns Real | Type check on return value | Wrong return type |
| Consistency | change_stat matches compute differences | Incorrect change_stat logic |

### Example Output

```julia
using ERGM, ERGMUserterms, Network

term = ExampleTerm()
net = Network{Int}(; n=10, directed=true)
for (i,j) in [(1,2),(2,3),(3,1),(1,4),(4,5)]
    add_edge!(net, i, j)
end

validate_term(term, net)
# [ Info: ✓ name() returns: example
# [ Info: ✓ compute() returns: 30.0 (Float64)
# [ Info: Testing change_stat() with 10 random edges...
# [ Info: ✓ change_stat() returns valid values
# [ Info: Running consistency check...
# [ Info: ✓ change_stat() consistent with compute()
# true
```

### Silent Validation

For programmatic use, disable verbose output:

```julia
if !validate_term(term, net; verbose=false)
    error("Term validation failed!")
end
```

### Early Termination

Validation stops at the first failure:

```julia
# If name() fails, compute() and change_stat() are not tested
# If compute() fails, change_stat() consistency is not tested
```

This prevents cascading errors and gives clear failure diagnostics.

## change_stat_check

Focused check that `change_stat()` correctly predicts differences in `compute()`:

```julia
consistent = change_stat_check(term, net; n_tests=10, verbose=true, tol=1e-10)
```

### How It Works

For each test:

1. Compute `stat_before = compute(term, net)`
2. Get `predicted = change_stat(term, net, i, j)`
3. Toggle edge `(i,j)` on a copy of the network
4. Compute `stat_after = compute(term, toggled_net)`
5. Check `|predicted - (stat_after - stat_before)| < tol`

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_tests` | Number of random edges to test | 10 |
| `verbose` | Print details of any failures | `true` |
| `tol` | Numerical tolerance | `1e-10` |

### Diagnosing Failures

When inconsistencies are found, verbose output shows the details:

```julia
change_stat_check(broken_term, net; verbose=true)
# [ Warning: Inconsistency at edge (3,7): predicted=2.0, actual=1.0
# false
```

Common causes of inconsistency:

| Cause | Symptom | Fix |
|-------|---------|-----|
| Wrong sign | predicted = -actual | Check `has_edge` direction in change_stat |
| Missing contribution | predicted < actual | Account for all affected edges/triangles |
| Double counting | predicted > actual | Avoid counting the same contribution twice |
| Off-by-one | close but not exact | Check edge iteration bounds |

### Increasing Test Coverage

```julia
# Quick check (default)
change_stat_check(term, net; n_tests=10)

# Thorough check
change_stat_check(term, net; n_tests=100)

# Very thorough
change_stat_check(term, net; n_tests=1000)
```

## consistency_check

A more thorough consistency check with an optional exhaustive mode:

```julia
# Random sampling (default)
consistent = consistency_check(term, net)

# Exhaustive check of ALL possible edges
consistent = consistency_check(term, net; exhaustive=true)
```

### Random vs Exhaustive

| Mode | Edges Checked | Speed | Completeness |
|------|---------------|-------|-------------|
| Random | min(100, n*n) random pairs | Fast | Probabilistic |
| Exhaustive | All n*(n-1) possible edges | Slow | Complete |

### When to Use Each Mode

```julia
# During development: use random for quick iteration
consistency_check(term, net)

# Before release: use exhaustive on small networks
small_net = Network{Int}(; n=10, directed=true)
# ... add edges ...
consistency_check(term, small_net; exhaustive=true)
```

### Numerical Tolerance

Both checking functions accept a tolerance parameter:

```julia
# Default tolerance (suitable for exact computations)
consistency_check(term, net; tol=1e-10)

# Relaxed tolerance (for floating-point intensive terms)
consistency_check(term, net; tol=1e-6)
```

## test_term

Comprehensive test suite that creates test networks and runs all validation checks:

```julia
passed = test_term(term; n_vertices=20, density=0.1, n_tests=100)
```

### What It Does

1. Creates a random directed network with specified size and density
2. Runs full `validate_term()` on the random network
3. Tests on an **empty network** (0 edges)
4. Tests on a **complete network** (all possible edges)
5. Reports results

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_vertices` | Size of random test network | 20 |
| `density` | Edge density (0.0 to 1.0) | 0.1 |
| `n_tests` | Random edge tests for consistency | 100 |

### Example Output

```text
Testing term: example
==================================================
[ Info: ✓ name() returns: example
[ Info: ✓ compute() returns: 232.0 (Float64)
[ Info: Testing change_stat() with 10 random edges...
[ Info: ✓ change_stat() returns valid values
[ Info: Running consistency check...
[ Info: ✓ change_stat() consistent with compute()

Additional tests:
✓ Works on empty network: 0.0
✓ Works on complete network: 4200.0
==================================================
All tests PASSED
```

### Testing Strategy

```julia
# Step 1: Quick validation during development
validate_term(term, net)

# Step 2: Comprehensive test when implementation is complete
test_term(term)

# Step 3: Test with different network sizes
test_term(term; n_vertices=5, density=0.3)   # Small dense
test_term(term; n_vertices=50, density=0.05)  # Large sparse
test_term(term; n_vertices=10, density=0.5)   # Medium moderate

# Step 4: Exhaustive consistency on small network
small_net = Network{Int}(; n=8, directed=true)
for _ in 1:15
    i, j = rand(1:8), rand(1:8)
    i != j && add_edge!(small_net, i, j)
end
consistency_check(term, small_net; exhaustive=true)
```

## Debugging Failed Validation

### Step 1: Identify the Failing Check

```julia
validate_term(term, net; verbose=true)
# Watch for ✗ or Warning messages
```

### Step 2: Test Methods Individually

```julia
# Test name
println("name: ", name(term))

# Test compute on a known network
net = Network{Int}(; n=3, directed=true)
add_edge!(net, 1, 2)
println("compute (1 edge): ", compute(term, net))

# Test change_stat
println("change_stat(1,2): ", change_stat(term, net, 1, 2))  # Existing edge
println("change_stat(2,3): ", change_stat(term, net, 2, 3))  # Non-existing edge
```

### Step 3: Manual Consistency Check

```julia
# Verify the fundamental relationship manually
net = Network{Int}(; n=4, directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)

before = compute(term, net)
delta = change_stat(term, net, 1, 3)

# Toggle edge and compute after
add_edge!(net, 1, 3)
after = compute(term, net)

println("Before: $before")
println("After: $after")
println("Predicted delta: $delta")
println("Actual delta: $(after - before)")
println("Match: $(abs(delta - (after - before)) < 1e-10)")
```

### Step 4: Check Edge Cases

```julia
# Empty network
empty_net = Network{Int}(; n=5)
println("Empty: ", compute(term, empty_net))

# Single edge
single_net = Network{Int}(; n=5, directed=true)
add_edge!(single_net, 1, 2)
println("Single edge: ", compute(term, single_net))

# Self-referential (if relevant)
println("change_stat(1,1): should handle gracefully")
```

## Best Practices

1. **Validate early and often**: Run `validate_term()` after each change to `compute()` or `change_stat()`
2. **Start with small networks**: Use 3-5 node networks where you can verify results by hand
3. **Test both edge states**: Verify change_stat for both existing and non-existing edges
4. **Use exhaustive checks**: Run `consistency_check(exhaustive=true)` on small networks before deployment
5. **Check numerical precision**: Use appropriate tolerance for floating-point computations
6. **Test multiple densities**: Validate on sparse, moderate, and dense networks
7. **Document expected values**: For small test cases, compute expected statistics by hand
8. **Automate testing**: Include validation in your test suite
