# Term Interface

This guide covers the ERGM term interface in detail, explaining how terms work and what your custom implementation must satisfy.

## The Three Required Methods

Every ERGM term must implement exactly three methods:

<!-- skip-check -->
```julia
name(term) -> String
compute(term, net) -> Float64
change_stat(term, net, i, j) -> Float64
```

These methods form a contract that ERGM.jl relies on for correct model estimation.

The examples below assume a term struct like this has been defined:

```julia
using ERGM, ERGMUserterms, Network
import ERGMUserterms: name, compute, change_stat

struct MyTerm <: AbstractUserTerm
    param::Float64
end

t = MyTerm(2.0)
net = network(10; directed=true)
```

### name()

Returns a human-readable string identifier for the term:

```julia
name(t::MyTerm) = "myterm"

# With parameters
name(t::MyTerm) = "myterm.$(t.param)"
```

#### Requirements

| Requirement | Description |
|-------------|-------------|
| Non-empty | Must return a non-empty string |
| Descriptive | Should identify the term and its configuration |
| Unique | Different configurations should produce different names |
| Deterministic | Must return the same value every time |

#### Naming Conventions

```julia
# Simple term
name(::Edges) = "edges"

# Parameterized term
name(t::GWDegree) = "gwdegree.$(t.decay)"

# Attribute-based term
name(t::NodeMatch) = "nodematch.$(t.attr)"

# Multiple parameters
name(t::InteractionTerm) = "interact.$(t.attr1).$(t.attr2)"
```

### compute()

Calculates the full network statistic value:

```julia
function compute(term::MyTerm, net)
    # Iterate over edges, vertices, etc.
    # Return the statistic value as a Float64
end
```

#### Requirements

| Requirement | Description |
|-------------|-------------|
| Returns Real | Must return a `Real` (typically `Float64`) |
| Deterministic | Same network must always produce the same value |
| No side effects | Must not modify the network |
| Handles edge cases | Should work on empty and complete networks |

#### Common Patterns

```julia
# Counting edges with a property
function compute(::MyEdgeTerm, net)
    count = 0.0
    for e in edges(net)
        # Check some condition
        count += 1.0
    end
    return count
end

# Summing over vertices
function compute(::MyNodeTerm, net)
    total = 0.0
    for v in vertices(net)
        total += length(outneighbors(net, v))^2
    end
    return total
end

# Using attributes. Note get_vertex_attribute returns an *empty Dict*
# (never `nothing`) when the attribute is absent — decide explicitly what
# missing values mean for your term (here they fall back to 0.0)
function compute(t::MyAttrTerm, net)
    attrs = get_vertex_attribute(net, t.attr)
    total = 0.0
    for e in edges(net)
        total += get(attrs, src(e), 0.0) * get(attrs, dst(e), 0.0)
    end
    return total
end
```

### change_stat()

Calculates the **add-direction** change statistic for dyad `(i,j)`:

```julia
function change_stat(term::MyTerm, net, i::Int, j::Int)
    # Return compute(net with edge i→j) - compute(net without edge i→j),
    # holding every other dyad fixed. The value must NOT depend on whether
    # the edge currently exists.
end
```

#### The Fundamental Relationship

```text
change_stat(term, net, i, j) == compute(term, net⁺ij) - compute(term, net⁻ij)
```

Where `net⁺ij`/`net⁻ij` are the network with edge `(i,j)` forced present/absent and all other dyads unchanged.

This relationship **must hold exactly**, for both current states of the dyad, for correct ERGM estimation: the MPLE design matrix uses the value directly and the Metropolis–Hastings sampler negates it for removal proposals. ERGMUserterms.jl validates both the value and its state-independence with [`change_stat_check`](@ref).

#### Requirements

| Requirement | Description |
|-------------|-------------|
| Returns Real | Must return a `Real` (typically `Float64`) |
| Consistent with compute | Must satisfy the fundamental relationship |
| State-independent | Same value whether the edge currently exists or not |
| Efficient | Should be O(degree), not O(edges) |

#### Sign Convention

Always return the *add-direction* value — do **not** flip the sign based on
the dyad's current state (the estimation machinery handles removals itself):

```julia
function change_stat(::MyTerm, net, i::Int, j::Int)
    # The contribution edge (i,j) makes when present. If the statistic
    # depends on degrees or adjacency, evaluate them with the dyad's own
    # edge masked out (the baseline state without i→j).
    return compute_local_contribution(net, i, j)
end
```

!!! warning "Common bug"
    The idiom `has_edge(net, i, j) ? -value : value` (toggle-direction) is
    **wrong** for ERGM.jl: it double-negates removal proposals in the MH
    sampler and sign-flips the MPLE design rows of observed edges.

#### Efficiency Guidelines

`change_stat()` is called thousands of times during MCMC simulation. It must be fast:

```julia
# GOOD: O(degree) - only examine neighbors of i and j
function change_stat(::TriangleTerm, net, i::Int, j::Int)
    shared = 0
    for k in outneighbors(net, i)
        k == j && continue
        has_edge(net, j, k) && (shared += 1)
    end
    return Float64(shared)
end

# BAD: O(n^2) - examines all pairs
function change_stat(::TriangleTerm, net, i::Int, j::Int)
    before = compute(TriangleTerm(), net)  # O(n^2)!
    # ... toggle and compute again
end
```

## Optional Trait: `is_dyad_dependent`

ERGM.jl classifies every term as dyad-dependent or dyad-independent with
`ERGM.is_dyad_dependent(term)`. The classification matters in two places:

- **MPLE honesty caveat**: pseudo-likelihood fits containing any
  dyad-dependent term print a standard-error warning in `show`.
- **MCMLE log-likelihood**: the bridge-sampling reference distribution
  zeroes the coefficients of dyad-dependent terms, so a misclassified term
  corrupts the reported log-likelihood/AIC/BIC.

For unknown term types the fallback is `true` — the conservative answer. If
your term's change statistic depends **only** on exogenous covariates of
dyad `(i,j)` (never on the state of other dyads), declare it:

<!-- skip-check -->
```julia
ERGM.is_dyad_dependent(::MyCovariateTerm) = false
```

!!! warning "Subtyping `NodalTerm`/`DyadicTerm` implies dyad-independence"
    ERGM.jl defines `is_dyad_dependent(::NodalTerm) = false` and
    `is_dyad_dependent(::DyadicTerm) = false`. Subtype those only for
    genuinely covariate-only terms; a term whose change statistic reads
    other dyads (degrees, shared partners, reciprocity) must not use them
    — keep `AbstractUserTerm` (fallback `true`) or add an explicit
    `is_dyad_dependent` method returning `true`.

## Attribute Validation Happens Only for Built-in Terms

`ERGMModel` construction validates ERGM.jl's *own* attribute-based terms
(`NodeCov`, `NodeMatch`, ...) against the network — a missing vertex
attribute throws an `ArgumentError` before any fitting — and snapshots
their attributes into typed vectors. User-defined terms pass through this
machinery **unchanged**: they are neither validated nor snapshotted. That
means:

- A user term reading a misspelled attribute will *not* be caught at model
  construction; decide explicitly how your `compute`/`change_stat` treat a
  missing attribute (error loudly, or document a default) and cover it with
  [`validate_term`](@ref)/[`change_stat_check`](@ref).
- For hot loops, do the typed-snapshot optimization yourself: read
  attributes once into a typed container in your term's constructor (or use
  Network.jl's typed accessors `vertex_attribute_vector(net, attr, V)` /
  `get_edge_attribute(net, attr, V)`) instead of hitting the untyped
  attribute Dicts in every `change_stat` call.

## Type Hierarchy

```text
AbstractERGMTerm (from ERGM.jl)
├── StructuralTerm
├── NodalTerm
├── DyadicTerm
└── AbstractUserTerm (from ERGMUserterms.jl)
```

### AbstractUserTerm

The base type for all user-defined terms:

<!-- skip-check -->
```julia
abstract type AbstractUserTerm <: AbstractERGMTerm end
```

Your terms should subtype `AbstractUserTerm`:

<!-- skip-check -->
```julia
struct MyTerm <: AbstractUserTerm
    # fields
end
```

This provides:

- Compatibility with all ERGM.jl functions
- Access to ERGMUserterms.jl validation tools
- Integration with the `@ergm_term` macro

### Choosing the Right Supertype

| Supertype | When to Use |
|-----------|-------------|
| `AbstractUserTerm` | Default choice for custom terms (`is_dyad_dependent` falls back to `true`) |
| `StructuralTerm` | Pure structural terms (no attributes) |
| `NodalTerm` | Covariate-only vertex-attribute terms (**implies dyad-independence**) |
| `DyadicTerm` | Covariate-only edge/dyadic-attribute terms (**implies dyad-independence**) |

For most custom terms, `AbstractUserTerm` is the right choice. Remember
that `NodalTerm`/`DyadicTerm` carry the `is_dyad_dependent = false` trait
(see above) — never use them for terms whose change statistic reads other
dyads.

## Struct Design

### Simple Terms (No Parameters)

```julia
struct MySimpleTerm <: AbstractUserTerm end
```

### Parameterized Terms

```julia
struct MyParamTerm <: AbstractUserTerm
    decay::Float64
    normalize::Bool
end

# Constructor with defaults
MyParamTerm(; decay=0.5, normalize=true) = MyParamTerm(decay, normalize)
```

### Generic Terms

```julia
struct MyGenericTerm{T} <: AbstractUserTerm
    param::T
    attr::Symbol
end

MyGenericTerm(param; attr=:none) = MyGenericTerm(param, attr)
```

### Attribute-Based Terms

```julia
struct MyAttrTerm <: AbstractUserTerm
    attr::Symbol
end

struct MyDyadTerm <: AbstractUserTerm
    covariate::Matrix{Float64}
end
```

## Working with Networks

### Accessing Network Structure

```julia
# Vertices
nv(net)              # Number of vertices
vertices(net)        # Iterator over vertex IDs

# Edges
ne(net)              # Number of edges
edges(net)           # Iterator over edges
has_edge(net, i, j)  # Check if edge exists

# Neighbors
outneighbors(net, v) # Outgoing neighbors
inneighbors(net, v)  # Incoming neighbors

# Edge endpoints
src(e)               # Source of edge
dst(e)               # Destination of edge
```

### Accessing Attributes

```julia
# Vertex attributes
attrs = get_vertex_attribute(net, :name)
# Returns Dict{Int, Any} or nothing

# Edge attributes
weights = get_edge_attribute(net, :weight)
# Returns Dict{Tuple{Int,Int}, Any} or nothing

# Network attributes
val = get_network_attribute(net, :name)
```

### Modifying Networks (Only in Tests)

```julia
# Add/remove edges
add_edge!(net, i, j)
rem_edge!(net, i, j)

# NEVER modify the network inside compute() or change_stat()
```

## Edge Cases to Handle

### Empty Networks

```julia
function compute(::MyTerm, net)
    ne(net) == 0 && return 0.0  # Handle empty network
    # ... normal computation
end
```

### Complete Networks

Your term should handle networks where all possible edges exist.

### Single-Vertex Networks

```julia
function compute(::MyTerm, net)
    nv(net) < 2 && return 0.0  # Need at least 2 vertices
    # ...
end
```

### Missing Attributes

```julia
function compute(t::MyAttrTerm, net)
    attrs = get_vertex_attribute(net, t.attr)
    isnothing(attrs) && return 0.0  # Attribute not set
    # ...
end
```

## Best Practices

1. **Implement compute() first**: Get the full statistic correct before optimizing change_stat()
2. **Use change_stat_check()**: Always verify consistency before using a term
3. **Avoid global state**: Terms should be pure functions of the network
4. **Keep structs immutable**: Don't use mutable state in term structs
5. **Return Float64**: Even for integer-valued statistics, return `Float64` for type stability
6. **Handle directed vs undirected**: Test your term on both network types
