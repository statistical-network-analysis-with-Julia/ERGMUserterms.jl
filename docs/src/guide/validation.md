# Validation and Testing

ERGMUserterms.jl provides a comprehensive validation framework to ensure custom terms are correctly implemented before use in ERGM estimation. Incorrect terms lead to biased estimates, so validation is critical.

## Overview

The validation process has three levels:

1. **Quick validation**: `validate_term()` checks all methods work and are consistent
2. **Consistency checking**: `change_stat_check()` and `consistency_check()` verify the fundamental relationship
3. **Comprehensive testing**: `test_term()` runs all checks on multiple network types

Every level draws its random dyads from one `rng::AbstractRNG` keyword
(default `Random.default_rng()`) and nothing else, so two calls with equal
rngs visit the same dyads and return the same verdict at any thread count —
and every failure message ends with a literal you can pass back to replay
it (see [Reproducing a failure](@ref)).

## validate_term

The primary validation function that checks all aspects of a term implementation:

<!-- skip-check -->
```julia
valid = validate_term(term, net; verbose=true, traits=true, n_tests=10,
                      rng=Random.default_rng())
```

### What It Checks

| Check | Description | Failure Means |
|-------|-------------|---------------|
| `name()` exists | Returns a non-empty string | Method not implemented |
| `name()` non-empty | String length > 0 | Empty name will cause display issues |
| `compute()` works | Returns without error | Bug in compute implementation |
| `compute()` returns Real | Type check on return value | Wrong return type |
| `change_stat()` works | Returns without error for `n_tests` random dyads | Bug in change_stat implementation |
| `change_stat()` returns Real | Type check on return value | Wrong return type |
| Consistency | change_stat matches compute differences on `n_tests` random dyads | Incorrect change_stat logic |
| Trait declarations | [`validate_traits`](@ref) (unless `traits=false`) | A declaration ERGM.jl would act on is wrong |

A check that could not run is never reported as ✓. In particular, an
attribute that is *constant* on `net` cannot be perturbed, so whether the
term reads it is undecidable there: `validate_traits` warns (`attribute :a
is constant on this network, so whether the term reads it is undecidable;
validate on a network where it varies`) and drops the ✓ line, and for an
attribute the term **declares** the run fails — validate on a network where
the attribute varies.

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `verbose` | Log each check and each failure | `true` |
| `traits` | Also run `validate_traits` | `true` |
| `n_tests` | Random dyads for the `change_stat` checks (at least 1; capped at the number of dyads) | 10 |
| `rng` | Source of every random draw; pass a `Xoshiro(seed)` for a reproducible run | `Random.default_rng()` |

`net` must be a network ERGM.jl can fit — at least 2 vertices, **one-mode**,
no self-loop — or it is refused with an `ArgumentError` before any check
runs; see [Networks the harness refuses](@ref). If `net`
has **masked dyads** and the term does not declare
`Networks.supports_missing`, the checks run on their face value (that is
what the term computes; ERGM.jl applies the `missing=` policy at estimation
time) and, with `verbose=true`, one `[ Info: net has k masked dyads; the term
declares supports_missing = false, so they are validated at their face value
…]` line says so — the same line appears on `change_stat_check` and on
`consistency_check(…; verbose=true)`.

### Worked Example

```julia
using ERGM, ERGMUserterms, Networks, Random

term = ExampleTerm()
net = network(10; directed=true)
for (i,j) in [(1,2),(2,3),(3,1),(1,4),(4,5)]
    add_edge!(net, i, j)
end

valid = validate_term(term, net; rng=Xoshiro(1))
# [ Info: ✓ name() returns: example
# [ Info: ✓ compute() returns: 26.0 (Float64)
# [ Info: Testing change_stat() with 10 random dyads...
# [ Info: ✓ change_stat() returns valid values
# [ Info: Running consistency check...
# [ Info: ✓ change_stat() consistent with compute()
# [ Info: ✓ declared attributes match the ones the term reads
# [ Info: ✓ is_dyad_dependent = false holds under toggling of other dyads
# [ Info: ✓ accepted by ERGMModel construction
# true
@assert valid
```

### Silent Validation

For programmatic use, run silently and check the returned flag — with
`verbose=false` nothing is logged, but the flag is still the full verdict:

```julia
if !validate_term(term, net; verbose=false, rng=Xoshiro(1))
    error("Term validation failed!")
end
```

### Early Termination

The interface checks all run, whatever fails: `name()`, `compute()` and the
`change_stat()` type check on `n_tests` random dyads each record their own
failure (with `verbose=true`, one warning each), so one run shows every
broken method. The two expensive stages are gated on everything before them
having passed:

- the **consistency check** (`change_stat` against brute-force `compute`
  differences) runs only if `name`, `compute` and the `change_stat` type
  check all passed — a `change_stat` that errors or returns a non-`Real`
  cannot be compared, and a `compute` that fails has nothing to compare it
  to;
- the **trait checks** ([`validate_traits`](@ref), unless `traits=false`)
  run only if the whole interface passed — the declarations of a term whose
  statistic is wrong are not worth exercising.

`verbose=false` changes only the logging; the returned flag is the same.

## change_stat_check

Focused check that `change_stat()` correctly predicts differences in `compute()`:

```julia
consistent = change_stat_check(term, net; n_tests=10, verbose=true, tol=1e-10,
                               rng=Random.default_rng())
@assert consistent
```

### How It Works

For each of `n_tests` random dyads `(i,j)` drawn from `rng` (never `i == j`;
on a copy of the network):

1. Get `predicted = change_stat(term, net, i, j)`
2. Compute `expected = compute(term, net⁺ij) − compute(term, net⁻ij)` by
   forcing the edge present/absent (the dyad's edge attributes are
   snapshotted and restored across the toggles)
3. Check `|predicted − expected| < tol`
4. **State-independence**: toggle the dyad and call
   `change_stat(term, net, i, j)` again — the value must not change
   (toggle-direction terms fail here)

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_tests` | Number of random dyads to test (at least 1: a run that checks no dyad returns no verdict) | 10 |
| `verbose` | Print details of any failures | `true` |
| `tol` | Numerical tolerance | `1e-10` |
| `rng` | Source of every random draw; pass a `Xoshiro(seed)` for a reproducible run | `Random.default_rng()` |

Every failure names the dyad and ends with the seed line
`reproduce with rng=Xoshiro(0x…, 0x…, 0x…, 0x…) passed to change_stat_check`.

### Diagnosing Failures

When inconsistencies are found, verbose output shows the details:

<!-- skip-check -->
```julia
change_stat_check(broken_term, net; verbose=true)
# [ Warning: Inconsistency at dyad (3,7): change_stat=2.0, brute-force=1.0;
#            reproduce with rng=Xoshiro(0x…, 0x…, 0x…, 0x…) passed to change_stat_check
# false
```

Toggle-direction terms (the pre-0.2 convention) fail the state-independence
step with a dedicated message:

<!-- skip-check -->
```julia
change_stat_check(old_convention_term, net; verbose=true)
# [ Warning: State-dependent change_stat at dyad (2,5): -1.0 with edge state
#            as-is vs 1.0 after toggling. change_stat must return the
#            add-direction change regardless of whether the edge exists;
#            reproduce with rng=Xoshiro(0x…, 0x…, 0x…, 0x…) passed to change_stat_check
# false
```

### Reproducing a failure

The `rng=Xoshiro(…)` literal at the end of every failure message is the
state of the harness's rng when the call started. The harness draws its
dyads straight from `rng`, so passing that literal back to the entry point
named in the message replays the identical dyad sequence — the failing dyad
included:

<!-- skip-check -->
```julia
change_stat_check(broken_term, net; verbose=true,
                  rng=Xoshiro(0x62bea62705299f9d, 0xb2196eb285598f6a,
                              0x813136b07248b437, 0x783171c16f41f41d))
```

The literal is exact for `Xoshiro` and for the default task-local rng (which
copies to one). For a `MersenneTwister` no short literal exists, so its
failures carry no hint — re-seed the generator you passed. `test_term`
prints its literal in the report header, since the failure may sit inside
the random network it generated.

Common causes of inconsistency:

| Cause | Symptom | Fix |
|-------|---------|-----|
| Toggle-direction sign flip | State-dependence warning; sign follows `has_edge` | Delete the `has_edge(net, i, j) ? -delta : delta` idiom; always return the add-direction value |
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
# Random sampling (default): distinct random dyads drawn from `rng`
consistent = consistency_check(term, net; rng=Xoshiro(1))
@assert consistent

# Exhaustive check of ALL possible edges (no randomness involved)
consistent = consistency_check(term, net; exhaustive=true)
@assert consistent
```

### Random vs Exhaustive

| Mode | Edges Checked | Speed | Completeness |
|------|---------------|-------|-------------|
| Random | up to min(100, number of dyads) distinct random dyads (never `i == j`; `(min, max)` on an undirected network) | Fast | Probabilistic |
| Exhaustive | Every dyad once: all n(n−1) ordered pairs on a directed network, all n(n−1)/2 unordered pairs on an undirected one | Slow | Complete |

`consistency_check` is silent by default; `verbose=true` reports the first
failing dyad with its replay literal. Both modes stop at the first failure.

### When to Use Each Mode

```julia
# During development: use random for quick iteration
consistency_check(term, net)

# Before release: use exhaustive on small networks
small_net = network(10; directed=true)
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
passed = test_term(term; n_vertices=20, density=0.1, n_tests=100,
                   rng=Xoshiro(1))
@assert passed
```

The generated networks are **directed and attribute-free unless you say
otherwise**. A term that declares `ERGM.required_vertex_attributes` needs
`vertex_attributes=`; one that declares `ERGM.requires_undirected` needs
`directed=false`. Without them `validate_term` reports the mismatch (`term
declares required vertex attribute :group, which the validation network does
not have` / `… validated on a directed network; validate it on an undirected
one`) and the run fails — correctly, but not because the term is wrong:

```julia
# The package template's term declares :group and requires_directed
using ERGM, ERGMUserterms, Random
include(joinpath(pkgdir(ERGMUserterms), "examples", "MyTermPackage", "src", "MyTermPackage.jl"))
using .MyTermPackage: ReciprocatedHomophily
@assert test_term(ReciprocatedHomophily(:group); n_vertices=12, n_tests=20, rng=Xoshiro(1),
                  vertex_attributes=Dict(:group => (v -> isodd(v) ? "a" : "b")))

# A term defined on undirected networks only is tested on undirected ones
@assert test_term(TemplateTerm(2.0); n_vertices=12, n_tests=20, directed=false, rng=Xoshiro(1))
```

### What It Does

1. Creates a random network with specified size and density (edges drawn
   from `rng`; directed unless `directed=false`) and sets the
   `vertex_attributes` on it
2. Runs full `validate_term()` on the random network with `n_tests` random
   dyads and the same `rng`
3. Tests on an **empty network** (0 edges) of the same kind
4. Tests on a **complete network** (all possible edges) of the same kind
5. Reports results

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_vertices` | Size of random test network (at least 2 — an `ArgumentError` otherwise) | 20 |
| `density` | Edge density, a proportion of the dyads (`0.0` to `1.0`; anything else is an `ArgumentError`) | 0.1 |
| `n_tests` | Random dyads for the `change_stat` checks (at least 1; passed to `validate_term`) | 100 |
| `directed` | Directed (`true`) or undirected (`false`) generated networks | `true` |
| `vertex_attributes` | `attr => spec` pairs set on every generated network; `spec` is a function `v -> value` or a vector with one value per vertex (length ≥ `n_vertices`) | `Dict{Symbol,Any}()` |
| `rng` | Source of every random draw (the network and the dyads) | `Random.default_rng()` |

`n_tests` is honoured: before 0.2.0 it was accepted and ignored (the network
was validated with a hard-wired 10 + 5 dyads); it is now forwarded to
`validate_term`, which caps it at the number of dyads in the generated
network. The report header prints the seed line
`(replay with rng=Xoshiro(0x…, 0x…, 0x…, 0x…))`, since a failure may sit
inside the random network `test_term` generated.

### Example Output

```text
Testing term: example
(replay with rng=Xoshiro(0x…, 0x…, 0x…, 0x…))
==================================================
[ Info: ✓ name() returns: example
[ Info: ✓ compute() returns: … (Float64)          # depends on the random network
[ Info: Testing change_stat() with 100 random dyads...
[ Info: ✓ change_stat() returns valid values
[ Info: Running consistency check...
[ Info: ✓ change_stat() consistent with compute()
[ Info: ✓ declared attributes match the ones the term reads
[ Info: ✓ is_dyad_dependent = false holds under toggling of other dyads
[ Info: ✓ accepted by ERGMModel construction

Additional tests:
✓ Works on empty network: 0.0
✓ Works on complete network: 990.0
==================================================
All tests PASSED
```

### Testing Strategy

```julia
# Step 1: Quick validation during development
validate_term(term, net)

# Step 2: Comprehensive test when implementation is complete
test_term(term)

# Step 3: Test with different network sizes (seeded, so a failure replays)
test_term(term; n_vertices=5, density=0.3, rng=Xoshiro(1))    # Small dense
test_term(term; n_vertices=50, density=0.05, rng=Xoshiro(2))  # Large sparse
test_term(term; n_vertices=10, density=0.5, rng=Xoshiro(3))   # Medium moderate

# Step 4: Exhaustive consistency on small network
rng = Xoshiro(4)
small_net = network(8; directed=true)
for _ in 1:15
    i, j = rand(rng, 1:8), rand(rng, 1:8)
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
net = network(3; directed=true)
add_edge!(net, 1, 2)
println("compute (1 edge): ", compute(term, net))

# Test change_stat
println("change_stat(1,2): ", change_stat(term, net, 1, 2))  # Existing edge
println("change_stat(2,3): ", change_stat(term, net, 2, 3))  # Non-existing edge
```

### Step 3: Manual Consistency Check

```julia
# Verify the fundamental relationship manually
net = network(4; directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)

before = compute(term, net)
delta = change_stat(term, net, 1, 3)

# Edge (1,3) is absent, so forcing it present gives the add-direction
# difference directly
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
empty_net = network(5)
println("Empty: ", compute(term, empty_net))

# Single edge
single_net = network(5; directed=true)
add_edge!(single_net, 1, 2)
println("Single edge: ", compute(term, single_net))

# Self-referential (if relevant)
println("change_stat(1,1): should handle gracefully")
```

## Networks the harness refuses

A network ERGM.jl cannot fit is refused **up front**, with an
`ArgumentError` carrying ERGM.jl's own reasoning, by every entry point that
takes a network — `validate_term`, `validate_traits`, `change_stat_check`,
`consistency_check` and `benchmark_term` — before any dyad is visited. It is
never validated dyad by dyad and then rejected at the `ERGMModel` step in
words that blame the term. Three cases:

- **Fewer than 2 vertices.** There is no dyad to check, so there is nothing
  a verdict could rest on. (`test_term` refuses `n_vertices < 2` and a
  `density` outside `[0, 1]` for the same reason, and every validator
  refuses `n_tests < 1` — a run that checks no dyad returns no verdict.)
- **Two-mode networks** (below).
- **A self-loop present.** ERGM.jl models the off-diagonal dyads only: a
  term's statistics would count a loop, but the pseudo-likelihood, `nobs`,
  the MH proposal and every simulation never touch the diagonal, so
  `ERGMModel` refuses the network (R ergm warns "This network contains
  loops" here). The harness would otherwise have validated `compute`
  *counting* the loop, every check green:

```julia
using ERGM, ERGMUserterms, Networks, Random
nl = network(6; directed=true, loops=true)
add_edge!(nl, 1, 1); add_edge!(nl, 1, 2)         # a self-loop at vertex 1
refused = try
    validate_term(ExampleTerm(), nl; rng=Xoshiro(1)); false
catch e
    e isa ArgumentError && occursin("contains 1 self-loop (at vertex 1)", e.msg)
end
@assert refused
rem_edge!(nl, 1, 1)                              # the fix the message suggests
@assert validate_term(ExampleTerm(), nl; verbose=false, rng=Xoshiro(1))
```

(`loops=true` is a capability, not a defect: a loops-allowed network with no
loop present validates normally.)

### Two-mode networks

ERGM.jl fits **one-mode** networks only: `ERGMModel` refuses a two-mode
(bipartite-flagged) network — `network(n; bipartite=k)` or a
`BipartiteNetwork` — because statnet's bipartite terms (`b1degree`,
`b2degree`, `b1factor`, `b1nodematch`, …) and the two-mode proposal kernel
are not implemented. The harness refuses one too, up front and with the same
message:

```julia
using ERGM, ERGMUserterms, Networks, Random
bp = network(6; directed=false, bipartite=3)     # modes {1,2,3} and {4,5,6}
add_edge!(bp, 1, 4)
@assert !add_edge!(bp, 1, 2)                     # a within-mode pair cannot hold an edge
refused = try
    validate_term(ExampleTerm(), bp; rng=Xoshiro(1)); false
catch e
    e isa ArgumentError && occursin("one-mode networks only", e.msg)
end
@assert refused
```

The refusal is not cosmetic. The harness draws and enumerates dyads over
*every* off-diagonal pair; on a two-mode network the within-mode pairs are
structurally absent (`add_edge!` returns `false`), so the brute-force
reference `compute(present) − compute(absent)` would read 0 there while a
correct `change_stat` returns its value, and a **correct term would be
reported inconsistent** (`Inconsistency at dyad (1,3): change_stat=4.0,
brute-force=0.0`). Validate on a one-mode network (`network(n; directed=…)`
without `bipartite=`); `test_term` and `profile_term` build one-mode,
loop-free networks of their own (directed unless `directed=false`, carrying
the `vertex_attributes` you pass — both take the same two keywords).

## Validation against R

The harness compares your `change_stat` with a brute-force reference
(`compute` with the dyad present minus `compute` with it absent). That
reference is only as good as the toggle-and-recompute code behind it, so it
is itself pinned against statnet: `test/fixtures/userterms_examples.toml`,
generated by the checked-in `test/fixtures/r/userterms_examples.R` (ergm
4.12.0, R 4.6.1, with a `[provenance]` block so
`Networks.load_golden` accepts it), records `summary()` and the complete
`ergmMPLE(output="array")` change-statistic array — every ordered dyad — for
the statnet terms the bundled terms are mathematically identical to:

| ERGMUserterms.jl term | statnet term |
|:--|:--|
| `ExampleTerm()` | `nodecov("id")`, `id` = vertex index |
| `TemplateTerm(1.0)` | `edges` |
| `WeightedEdges()` | `edgecov(W)`, `W` = 2.5 on the weighted arcs and 1.0 (the `default`) elsewhere |
| `DyadCovTerm(cov)` | `edgecov(cov)`, asymmetric `cov[i,j] = 2i + j` |
| `InteractionTerm(:a, :b)` | `edgecov(M)`, `M[i,j] = a_i b_j + a_j b_i` |
| `ReciprocatedHomophily(:group)` (`examples/MyTermPackage`) | `mutual(same="group", diff=FALSE)` |

The network is an 8-vertex directed graph from a literal edge list (the
template's own test network plus four arcs, so every change statistic takes
more than one value) and its undirected projection, on which `W` and `cov`
are symmetrised to their `(min, max)` entry — the canonical key
`WeightedEdges`/`DyadCovTerm` read on an undirected network, which is how
that branch gets an R number too. The testset "Golden: bundled terms and
the harness match statnet" compares, at the fixture's `1e-9` (deterministic
integer/rational arithmetic on a fixed graph — no estimator, no Monte Carlo,
so any disagreement is a term bug), each term's `compute`, each term's
`change_stat` on every dyad (in both dyad orders when undirected), and
`ERGMUserterms._brute_change_stat` on every dyad. R's coefficient labels
(`nodecov.id`, `edgecov.W`, …) are recorded in the fixture but deliberately
*not* claimed: a user term keeps the name its author gives it.

Regenerate the fixture from the package root with
`Rscript test/fixtures/r/userterms_examples.R > test/fixtures/userterms_examples.toml`
(needs `ergm` installed in R).

## Best Practices

1. **Validate early and often**: Run `validate_term()` after each change to `compute()` or `change_stat()`
2. **Start with small networks**: Use 3-5 node networks where you can verify results by hand
3. **Test both edge states**: Verify change_stat for both existing and non-existing edges
4. **Use exhaustive checks**: Run `consistency_check(exhaustive=true)` on small networks before deployment
5. **Check numerical precision**: Use appropriate tolerance for floating-point computations
6. **Test multiple densities**: Validate on sparse, moderate, and dense networks
7. **Document expected values**: For small test cases, compute expected statistics by hand
8. **Automate testing**: Include validation in your test suite, with an explicit `rng=Xoshiro(seed)` so a red run is the same run every time
