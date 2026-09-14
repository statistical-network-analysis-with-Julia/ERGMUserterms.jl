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
using ERGM, ERGMUserterms, Networks, Random
using Graphs: src, dst                    # edge endpoints are Graphs.jl's
import ERGM: name, compute, change_stat   # the shared generics (ERGMUserterms re-exports the same ones)

struct MyTerm <: AbstractUserTerm
    param::Float64
end

t = MyTerm(2.0)
net = network(10; directed=true)
```

The `import` line is not optional: `name`, `compute` and `change_stat` are
ERGM.jl's generics (`compute` and `name` are Networks.jl's shared statistic
protocol), and ERGMUserterms.jl re-exports the very same functions —
`ERGMUserterms.compute === ERGM.compute === Networks.compute`. A method
added after a bare `using` defines a rival local function that ERGM.jl never
calls. `src`/`dst` are `Graphs.src`/`Graphs.dst`, which neither Networks.jl
nor ERGM.jl re-exports, so a term iterating `edges(net)` imports them.

### name()

Returns a human-readable string identifier for the term — `"myterm"` for a
term without parameters, and with the parameters spelled out otherwise:

```julia
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

<!-- skip-check -->
```julia
# Simple term
name(::Edges) = "edges"

# Parameterized term (R ergm's label, so coefficient tables line up with statnet)
name(t::GWDegree) = "gwdeg.fixed.$(t.decay)"

# Attribute-based term
name(t::NodeMatch) = "nodematch.$(t.attr)"

# Multiple parameters
name(t::InteractionTerm) = "interact.$(t.attr1).$(t.attr2)"
```

ERGM.jl also has a two-argument `name(term, net)` for labels that depend on
the network's direction (`name(GWESP(0.5), net)` is `"gwesp.fixed.0.5"` on an
undirected network and `"gwesp.OTP.fixed.0.5"` on a directed one, as in
statnet). It falls back to `name(term)`, so a user term needs only the
one-argument method; add the two-argument one only if your label really
depends on the network.

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

<!-- skip-check -->
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

# Using attributes. get_vertex_attribute returns an *empty Dict* (never
# `nothing`) when the attribute is absent. DECLARE the attribute
# (`ERGM.required_vertex_attributes(t::MyAttrTerm) = (t.attr,)`): then an
# ERGMModel refuses a network on which it is absent or set on only some
# vertices, and the `get(attrs, v, 0.0)` fallback below is reached only by
# raw compute calls — see "Missing Attribute Values (NA)" at the bottom
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

<!-- skip-check -->
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
    sampler and sign-flips the MPLE design rows of observed edges. It is,
    however, R ergm's own C convention — `C_CHANGESTAT_FN(…, Rboolean
    edgestate)` in `ergm_changestat.h`, with every statnet changestat
    written `CHANGE_STAT[0] += edgestate ? -1 : 1` — so it is exactly what
    a port of an `ergm.userterms` changestat arrives with. When porting,
    **drop the `edgestate` sign and return the add-direction value**; ERGM.jl
    applies the sign for removals itself.

#### Efficiency Guidelines

`change_stat()` is called thousands of times during MCMC simulation. It must be fast:

<!-- skip-check -->
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

## Declaring Your Term's Traits

Implementing `name`/`compute`/`change_stat` makes a term *computable*.
Declaring its **traits** makes it a first-class citizen: ERGM.jl's formula
validation and materialization act on the trait declarations, not on ERGM's
own term types, so a term from any package is checked exactly like a built-in
one at `ERGMModel` construction.

| trait | default | declare it when |
|:------|:--------|:----------------|
| `ERGM.required_vertex_attributes(t)` | `()` | the term reads a vertex attribute and is meaningless without it |
| `ERGM.required_edge_attributes(t)` | `()` | likewise for an edge attribute |
| `ERGM.requires_directed(t)` | `false` | the statistic is undefined on undirected networks (ERGM's `Mutual`, `OStar`, `GWODegree`, … declare it) |
| `ERGM.requires_undirected(t)` | `false` | ... or on directed ones (ERGM's `Degree`, `Kstar`, `GWDegree` declare it) |
| `ERGM.is_dyad_dependent(t)` | `true` | the change statistic never reads another dyad → declare `false` |
| `Networks.supports_missing(t)` | `false` | the statistic consults `is_missing_dyad` and ignores masked dyads (`ERGM.supports_missing` is the same function) |

<!-- skip-check -->
```julia
struct Homophily <: AbstractUserTerm
    attr::Symbol
end

ERGM.required_vertex_attributes(t::Homophily) = (t.attr,)
ERGM.is_dyad_dependent(::Homophily) = false
```

With that declaration, `ERGMModel(ERGMFormula([Edges(), Homophily(:groop)]), net)`
throws an `ArgumentError` naming the misspelled attribute. Without it, the term
would read an empty attribute Dict, evaluate to a constant zero, and produce a
meaningless coefficient — silently.

!!! warning "Declare an attribute only if its absence is an error"
    "Required" means *the model cannot be fitted without it*. A term that reads
    an attribute **if present** and falls back to a default — like
    [`WeightedEdges`](@ref), whose edges without a `:weight` count at
    `t.default` — genuinely does not require it, and declaring it would reject
    networks the term handles perfectly well.

[`validate_term`](@ref) exercises every declaration (see
[`validate_traits`](@ref)): it perturbs the network's attributes to check that
the ones you declared are the ones your term actually reads, checks that a
direction requirement is enforced by `ERGMModel`, toggles other dyads to check
a dyad-independence claim, and flips the face value of masked dyads to check a
`supports_missing` claim. The package template in
`examples/MyTermPackage/` is a complete term declaring all four.

### `supports_missing`: the same idiom as ERGM's own terms

`Networks.supports_missing` and `ERGM.supports_missing` are one function
(`ERGM.supports_missing === Networks.supports_missing`), so either spelling
extends the same generic. Declare `true` only if the statistic consults the
mask itself — the same `ObservedEdges` idiom ERGM.jl's own `supports_missing`
docstring uses:

```julia
struct ObservedEdges <: AbstractUserTerm end
name(::ObservedEdges) = "observed_edges"
function compute(::ObservedEdges, net)
    total = 0.0
    for e in edges(net)
        is_missing_dyad(net, src(e), dst(e)) && continue   # the mask, not the face value
        total += 1.0
    end
    return total
end
change_stat(::ObservedEdges, net, i::Int, j::Int) = is_missing_dyad(net, i, j) ? 0.0 : 1.0
Networks.supports_missing(::ObservedEdges) = true    # ≡ ERGM.supports_missing(::ObservedEdges) = true

masked = network(5; directed=true)
add_edge!(masked, 1, 2); add_edge!(masked, 2, 3)
set_missing_dyad!(masked, 1, 2)
compute(ObservedEdges(), masked)           # 1.0 — the masked tie does not count
@assert compute(ObservedEdges(), masked) == 1.0
@assert validate_traits(ObservedEdges(), masked; verbose=false, rng=Xoshiro(1))
```

`validate_traits` checks the claim by flipping the face value of a masked
dyad: the statistic of a term declaring `true` must not move. The default
`false` is the honest answer for a term that counts edges as stored — that is
what every ERGM.jl built-in declares, and the missing-data treatment then
lives in the estimator (MPLE drops masked dyads from the pseudo-likelihood).

### `is_dyad_dependent` in detail

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

## Attribute Validation and Snapshotting

`ERGMModel` construction validates every term — built-in or not — against the
network, using the attributes the term **declares** (see
[Declaring Your Term's Traits](@ref "Declaring Your Term's Traits") above): a
declared vertex or edge attribute that the network does not have throws an
`ArgumentError` before any fitting, and so does a declared *vertex* attribute
with a value on only some vertices (statnet refuses NA — see
[Missing Attribute Values (NA)](@ref) below). The check is
`ERGM._validate_formula`, public and documented in ERGM.jl. A term that
declares nothing is not validated, so an undeclared attribute-based term will
read an empty Dict and evaluate to zero — declare your attributes.

Typed *snapshotting* (materialization into dense `Vector{Float64}`/code
vectors) is still done only for ERGM.jl's own nodal terms. For hot loops, do
that optimization yourself: read attributes once into a typed container in your
term's constructor (or use Networks.jl's typed accessors
`vertex_attribute_vector(net, attr, V)` / `get_edge_attribute(net, attr, V)`).
Short of that, read per element — `get_vertex_attribute(net, attr, v)` /
`get_edge_attribute(net, attr, i, j)` return `nothing` when the value is absent
and allocate nothing — and assert the value's type (`Float64(x)::Float64`), as
the bundled `WeightedEdges`/`InteractionTerm` do; never `get(…)` from the
whole-Dict getters in `change_stat`, whose `Dict{…,Any}` makes the method
infer `Any` and box on every sampler step. (The typed snapshot is faster
still: a vector index instead of a hash lookup in every `change_stat` call.)

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

A supertype does not fix a term's direction requirement; a declaration does.
ERGM.jl's own `Degree`, `Kstar` and `GWDegree` are the built-in examples of
`requires_undirected = true` — exactly the terms R ergm refuses on a
directed network — and they refuse a directed network with one message
naming the directed variant to use instead: `ERGMModel` (and `compute` /
`change_stat` themselves) throw `term 'kstar2' is only defined for
undirected networks, but the network is directed. Use OStar(2) / IStar(2)
…` (`GWODegree`/`GWIDegree` for `GWDegree`, `ODegree`/`IDegree` for
`Degree`). A user term declaring `ERGM.requires_undirected(::MyTerm) = true`
gets the same refusal from `ERGMModel` — that is the built-in analogue of
the declaration, and [`validate_traits`](@ref) checks it is enforced.

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

<!-- skip-check -->
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

<!-- skip-check -->
```julia
# Vertex attributes
attrs = get_vertex_attribute(net, :name)
# Returns Dict{Int, Any}; an *empty* Dict (never nothing) when the
# attribute is absent

# Edge attributes
weights = get_edge_attribute(net, :weight)
# Returns Dict{Tuple{Int,Int}, Any}; an *empty* Dict (never nothing) when
# the attribute is absent. Keys are canonical: (i,j) for directed
# networks, (min,max) for undirected

# Network attributes
val = get_network_attribute(net, :name)
```

### Modifying Networks (Only in Tests)

<!-- skip-check -->
```julia
# Add/remove edges
add_edge!(net, i, j)
rem_edge!(net, i, j)

# NEVER modify the network inside compute() or change_stat()
```

## Edge Cases to Handle

### Empty Networks

<!-- skip-check -->
```julia
function compute(::MyTerm, net)
    ne(net) == 0 && return 0.0  # Handle empty network
    # ... normal computation
end
```

### Complete Networks

Your term should handle networks where all possible edges exist.

### Single-Vertex Networks

<!-- skip-check -->
```julia
function compute(::MyTerm, net)
    nv(net) < 2 && return 0.0  # Need at least 2 vertices
    # ...
end
```

(`compute` may be called on such a network directly; the harness itself
refuses one — with fewer than 2 vertices there is no dyad to check, so
`validate_term` and friends throw an `ArgumentError` rather than report a
verdict resting on nothing.)

### Missing Attribute Values (NA)

`get_vertex_attribute`/`get_edge_attribute` return an *empty* Dict (never
`nothing`) when the attribute is absent, and vertices may lack individual
values. What that means for your term is settled by the contract, not by the
term:

- A **declared** vertex attribute (`ERGM.required_vertex_attributes`) must be
  set on **every** vertex, or `ERGMModel` construction throws — the message
  is ERGM's own: `term 'attrsum.a' needs vertex attribute :a on every vertex,
  but 1 of 4 vertices have no value (vertices 4). statnet refuses NA
  attribute values; set a value for every vertex or drop the term.` The check
  is `ERGM._validate_formula` (public, documented in ERGM.jl), the same one
  that validates ERGM's built-in attribute terms.
- An **undeclared** attribute read is a [`validate_traits`](@ref) failure:
  on a network lacking the attribute the term would silently become an
  all-zero design column.
- A `get(attrs, v, default)` fallback in your code is therefore only ever
  reached by raw `compute`/`change_stat` calls outside a model. Keep it (it
  makes the term total on any network), but do not rely on it to define what
  a missing value means.

```julia
struct MyAttrSum <: AbstractUserTerm
    attr::Symbol
end
name(t::MyAttrSum) = "attrsum.$(t.attr)"
# One vertex's value, read per element (0 B, infers Float64); the zero-fill
# is reached only by raw compute calls — inside an ERGMModel the declared
# attribute is complete
function _attr_value(t::MyAttrSum, net, v::Int)
    x = get_vertex_attribute(net, t.attr, v)
    return x === nothing ? 0.0 : Float64(x)::Float64
end
function compute(t::MyAttrSum, net)
    total = 0.0
    for e in edges(net)
        total += _attr_value(t, net, Int(src(e))) + _attr_value(t, net, Int(dst(e)))
    end
    return total
end
change_stat(t::MyAttrSum, net, i::Int, j::Int) =
    _attr_value(t, net, i) + _attr_value(t, net, j)
ERGM.required_vertex_attributes(t::MyAttrSum) = (t.attr,)
ERGM.is_dyad_dependent(::MyAttrSum) = false

partial = network(4; directed=true)
add_edge!(partial, 1, 2)
set_vertex_attribute!(partial, :a, Dict(1 => 1.0, 2 => 2.0, 3 => 3.0))   # vertex 4: no value
compute(MyAttrSum(:a), partial)                       # 3.0 — the raw call zero-fills
err = try
    ERGMModel(ERGMFormula([MyAttrSum(:a)]), partial)  # the model refuses
    nothing
catch e
    e
end
@assert err isa ArgumentError
@assert occursin("on every vertex", err.msg)
# validate_traits reports the vertex without a value (quietly here) and fails
@assert !validate_traits(MyAttrSum(:a), partial; verbose=false, rng=Xoshiro(1))

set_vertex_attribute!(partial, :a, 4, 4.0)             # complete it ...
@assert validate_term(MyAttrSum(:a), partial; verbose=false, rng=Xoshiro(1))   # ... and it is valid
```

## Best Practices

1. **Implement compute() first**: Get the full statistic correct before optimizing change_stat()
2. **Use change_stat_check()**: Always verify consistency before using a term
3. **Avoid global state**: Terms should be pure functions of the network
4. **Keep structs immutable**: Don't use mutable state in term structs
5. **Return Float64**: Even for integer-valued statistics, return `Float64` for type stability
6. **Handle directed vs undirected**: Test your term on both network types
