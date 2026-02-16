# Term Interface

This guide covers the ERGM term interface in detail, explaining how terms work and what your custom implementation must satisfy.

## The Three Required Methods

Every ERGM term must implement exactly three methods:

```julia
name(term) -> String
compute(term, net) -> Float64
change_stat(term, net, i, j) -> Float64
```

These methods form a contract that ERGM.jl relies on for correct model estimation.

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

# Using attributes
function compute(t::MyAttrTerm, net)
    attrs = get_vertex_attribute(net, t.attr)
    isnothing(attrs) && return 0.0
    total = 0.0
    for e in edges(net)
        total += get(attrs, src(e), 0.0) * get(attrs, dst(e), 0.0)
    end
    return total
end
```

### change_stat()

Calculates the change in the statistic when edge `(i,j)` is toggled:

```julia
function change_stat(term::MyTerm, net, i::Int, j::Int)
    # If edge exists: return compute(net_without_ij) - compute(net)
    # If edge absent: return compute(net_with_ij) - compute(net)
end
```

#### The Fundamental Relationship

```text
change_stat(term, net, i, j) == compute(term, net') - compute(term, net)
```

Where `net'` is the network with edge `(i,j)` toggled (added if absent, removed if present).

This relationship **must hold exactly** for correct ERGM estimation. ERGMUserterms.jl validates this with [`change_stat_check`](@ref).

#### Requirements

| Requirement | Description |
|-------------|-------------|
| Returns Real | Must return a `Real` (typically `Float64`) |
| Consistent with compute | Must satisfy the fundamental relationship |
| Handles both directions | Must work whether edge exists or not |
| Efficient | Should be O(degree), not O(edges) |

#### Sign Convention

The sign depends on whether the edge currently exists:

```julia
function change_stat(::MyTerm, net, i::Int, j::Int)
    value = compute_local_contribution(net, i, j)
    return has_edge(net, i, j) ? -value : value
end
```

- **Edge absent → being added**: Return the positive contribution
- **Edge present → being removed**: Return the negative contribution

#### Efficiency Guidelines

`change_stat()` is called thousands of times during MCMC simulation. It must be fast:

```julia
# GOOD: O(degree) - only examine neighbors of i and j
function change_stat(::TriangleTerm, net, i::Int, j::Int)
    shared = 0
    for k in outneighbors(net, i)
        has_edge(net, j, k) && (shared += 1)
    end
    return has_edge(net, i, j) ? -Float64(shared) : Float64(shared)
end

# BAD: O(n^2) - examines all pairs
function change_stat(::TriangleTerm, net, i::Int, j::Int)
    before = compute(TriangleTerm(), net)  # O(n^2)!
    # ... toggle and compute again
end
```

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

```julia
abstract type AbstractUserTerm <: AbstractERGMTerm end
```

Your terms should subtype `AbstractUserTerm`:

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
| `AbstractUserTerm` | Default choice for custom terms |
| `StructuralTerm` | Pure structural terms (no attributes) |
| `NodalTerm` | Terms using vertex attributes |
| `DyadicTerm` | Terms using edge/dyadic attributes |

For most custom terms, `AbstractUserTerm` is the right choice.

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
