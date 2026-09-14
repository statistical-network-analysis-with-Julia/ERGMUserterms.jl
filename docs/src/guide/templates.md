# Templates and Examples

ERGMUserterms.jl provides five complete example terms that serve as templates for building your own. Each demonstrates a different pattern for custom term development.

## Template Overview

| Template | Pattern | Key Feature |
|----------|---------|-------------|
| [`ExampleTerm`](@ref) | Simple structural | No parameters, basic edge iteration |
| [`TemplateTerm`](@ref) | Parameterized | Generic type parameter, keyword constructor |
| [`WeightedEdges`](@ref) | Edge attributes | Accessing edge-level covariates |
| [`DyadCovTerm`](@ref) | Dyadic covariate | Using an external matrix |
| [`InteractionTerm`](@ref) | Node interactions | Combining multiple vertex attributes |

## Package Template: `examples/MyTermPackage`

The five terms above are *term* templates. For shipping terms in a package of
your own, the repository also carries a **package template** —
`examples/MyTermPackage/` — a copyable skeleton (`Project.toml`, module,
tests) whose term `ReciprocatedHomophily` declares the full set of ERGM.jl
traits:

<!-- skip-check -->
```julia
ERGM.required_vertex_attributes(t::ReciprocatedHomophily) = (t.attr,)
ERGM.requires_directed(::ReciprocatedHomophily)           = true
ERGM.is_dyad_dependent(::ReciprocatedHomophily)           = true
Networks.supports_missing(::ReciprocatedHomophily)        = true
```

and, on the strength of those declarations, is validated by
[`validate_term`](@ref) and accepted (or correctly *rejected* — on an
undirected network, or one lacking the attribute) by `ERGMModel` construction,
exactly like an ERGM.jl built-in. See
[Declaring Your Term's Traits](@ref "Declaring Your Term's Traits").

## ExampleTerm

The simplest possible custom term. Counts edges weighted by vertex ID sum.

### Implementation

The template implementations below re-define local copies of the bundled
terms so you can adapt them:

```julia
using ERGM, ERGMUserterms, Networks, Random
using Graphs: src, dst                     # edge endpoints are Graphs.jl's
import ERGM: name, compute, change_stat    # the shared generics

struct ExampleTerm <: AbstractUserTerm end

name(::ExampleTerm) = "example"

function compute(::ExampleTerm, net)
    total = 0.0
    for e in edges(net)
        total += src(e) + dst(e)
    end
    return total
end

function change_stat(::ExampleTerm, net, i::Int, j::Int)
    # Add-direction: adding edge (i,j) adds i + j, whatever the current state
    return Float64(i + j)
end

# Covariate-only: opt out of ERGM.jl's conservative dyad-dependent fallback
ERGM.is_dyad_dependent(::ExampleTerm) = false
```

### Usage

```julia
term = ExampleTerm()
net = network(5; directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 3, 4)

compute(term, net)  # (1+2) + (3+4) = 10.0
change_stat(term, net, 2, 3)  # +(2+3) = 5.0
change_stat(term, net, 1, 2)  # +(1+2) = 3.0 (state-independent add-direction)
@assert compute(term, net) == 10.0
@assert validate_term(term, net; verbose=false, rng=Xoshiro(1))
```

### When to Use This Pattern

- Structural terms with no parameters
- Terms that iterate over all edges
- Terms where the change statistic is a simple local computation

### Adapting This Template

```julia
# Example: Count edges where source ID > target ID
struct DescendingEdges <: AbstractUserTerm end

name(::DescendingEdges) = "descending"

function compute(::DescendingEdges, net)
    count = 0.0
    for e in edges(net)
        src(e) > dst(e) && (count += 1.0)
    end
    return count
end

function change_stat(::DescendingEdges, net, i::Int, j::Int)
    # Add-direction: 1 if this edge would contribute, else 0
    return i > j ? 1.0 : 0.0
end
ERGM.is_dyad_dependent(::DescendingEdges) = false

@assert validate_term(DescendingEdges(), net; verbose=false, rng=Xoshiro(1))
```

## TemplateTerm

A parameterized template with a generic type parameter. Copy and modify for terms that need configuration.

### Implementation

```julia
struct TemplateTerm{T} <: AbstractUserTerm
    param::T
    attr::Symbol

    TemplateTerm(param::T; attr::Symbol=:none) where T = new{T}(param, attr)
end

name(t::TemplateTerm) = "template.$(t.param)"

function compute(t::TemplateTerm, net)
    # Count edges multiplied by parameter
    return Float64(ne(net)) * t.param
end

function change_stat(t::TemplateTerm, net, i::Int, j::Int)
    # Add-direction: adding edge (i,j) increases the statistic by the
    # parameter (state-independent; removals are handled by the sampler)
    return Float64(t.param)
end

# Covariate-only term: opt out of ERGM.jl's conservative dyad-dependent
# fallback (omit this line for terms that read the state of other dyads)
ERGM.is_dyad_dependent(::TemplateTerm) = false
```

### Usage

```julia
# Create with different parameter types
term_float = TemplateTerm(2.5)
term_int = TemplateTerm(3)

net = network(5; directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)

compute(term_float, net)  # 2 * 2.5 = 5.0
change_stat(term_float, net, 3, 4)  # +2.5 (edge absent)
change_stat(term_float, net, 1, 2)  # +2.5 (edge present: same add-direction value)
@assert validate_term(term_float, net; verbose=false, rng=Xoshiro(1))
```

### When to Use This Pattern

- Terms with a scalar parameter
- Terms where the parameter scales the statistic linearly
- Prototype terms before implementing efficient change_stat()

### Adapting This Template

```julia
# Example: Weighted edge count with decay
struct DecayEdges <: AbstractUserTerm
    decay::Float64
    DecayEdges(; decay::Float64=0.5) = new(decay)
end

name(t::DecayEdges) = "decay_edges.$(t.decay)"

function compute(t::DecayEdges, net)
    total = 0.0
    for e in edges(net)
        # Weight by position in vertex ordering
        total += t.decay ^ abs(src(e) - dst(e))
    end
    return total
end

function change_stat(t::DecayEdges, net, i::Int, j::Int)
    return t.decay ^ abs(i - j)
end
ERGM.is_dyad_dependent(::DecayEdges) = false

@assert validate_term(DecayEdges(), net; verbose=false, rng=Xoshiro(1))
```

## WeightedEdges

Demonstrates accessing edge attributes from the network.

### Implementation

```julia
struct WeightedEdges <: AbstractUserTerm
    attr::Symbol
    default::Float64
    WeightedEdges(attr::Symbol=:weight; default::Float64=1.0) = new(attr, default)
end

name(t::WeightedEdges) = "weightedges.$(t.attr)"

# The weight edge (i,j) carries: its stored value, else the default a fresh
# edge would get. The per-edge getter canonicalises the key on undirected
# networks ((min,max)) and returns `nothing` for an edge — or a network —
# without the attribute; the `::Float64` assertion keeps the term type-stable
# and allocation-free (the attribute store is an untyped `Dict{…,Any}`, and a
# bare `get` from it would infer `Any` and box on every sampler step)
function _edge_weight(t::WeightedEdges, net, i::Int, j::Int)
    w = get_edge_attribute(net, t.attr, i, j)
    return w === nothing ? t.default : Float64(w)::Float64
end

function compute(t::WeightedEdges, net)
    # Sum over the network's edges with a default — this keeps compute
    # consistent with change_stat when the sampler adds unweighted edges
    total = 0.0
    for e in edges(net)
        total += _edge_weight(t, net, Int(src(e)), Int(dst(e)))
    end
    return total
end

function change_stat(t::WeightedEdges, net, i::Int, j::Int)
    # Add-direction: the weight edge (i,j) carries, whatever the dyad's state
    return _edge_weight(t, net, i, j)
end
ERGM.is_dyad_dependent(::WeightedEdges) = false
# Deliberately NO `ERGM.required_edge_attributes`: edges without the
# attribute count at the default, so the term is well defined without it
```

### Usage

```julia
net = network(5; directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)
set_edge_attribute!(net, :weight, 1, 2, 3.0)
set_edge_attribute!(net, :weight, 2, 3, 1.5)

term = WeightedEdges(:weight)
compute(term, net)  # 3.0 + 1.5 = 4.5
@assert compute(term, net) == 4.5
@assert validate_term(term, net; verbose=false, rng=Xoshiro(1))
```

### When to Use This Pattern

- Terms that use edge-level data (weights, types, timestamps)
- Terms with a configurable attribute name
- Terms that fall back to a default when attributes are missing

### Adapting This Template

```julia
# Example: Sum of squared edge weights
struct SquaredWeights <: AbstractUserTerm
    attr::Symbol
    SquaredWeights(attr::Symbol=:weight) = new(attr)
end

name(t::SquaredWeights) = "sqweights.$(t.attr)"

function _weight(t::SquaredWeights, net, i::Int, j::Int)
    w = get_edge_attribute(net, t.attr, i, j)
    return w === nothing ? 1.0 : Float64(w)::Float64
end

function compute(t::SquaredWeights, net)
    total = 0.0
    for e in edges(net)
        total += _weight(t, net, Int(src(e)), Int(dst(e)))^2
    end
    return total
end

function change_stat(t::SquaredWeights, net, i::Int, j::Int)
    return _weight(t, net, i, j)^2
end
ERGM.is_dyad_dependent(::SquaredWeights) = false

@assert validate_term(SquaredWeights(), net; verbose=false, rng=Xoshiro(1))
```

## DyadCovTerm

Demonstrates using an external matrix as a dyadic covariate.

### Implementation

```julia
struct DyadCovTerm <: AbstractUserTerm
    covariate::Matrix{Float64}
end

name(::DyadCovTerm) = "dyadcov"

function compute(t::DyadCovTerm, net)
    total = 0.0
    for e in edges(net)
        i, j = src(e), dst(e)
        if i <= size(t.covariate, 1) && j <= size(t.covariate, 2)
            total += t.covariate[i, j]
        end
    end
    return total
end

function change_stat(t::DyadCovTerm, net, i::Int, j::Int)
    # Undirected edges are stored canonically as (min, max); read the
    # covariate the same way compute() sees the edge
    if !is_directed(net)
        i, j = minmax(i, j)
    end
    i <= size(t.covariate, 1) && j <= size(t.covariate, 2) || return 0.0
    return t.covariate[i, j]
end
ERGM.is_dyad_dependent(::DyadCovTerm) = false
```

### Usage

```julia
# Distance matrix between actors
distance = [
    0.0  1.0  2.0  3.0
    1.0  0.0  1.5  2.5
    2.0  1.5  0.0  1.0
    3.0  2.5  1.0  0.0
]

net = network(4; directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 3, 4)

term = DyadCovTerm(distance)
compute(term, net)  # distance[1,2] + distance[3,4] = 1.0 + 1.0 = 2.0
change_stat(term, net, 1, 3)  # +distance[1,3] = 2.0 (adding edge)
@assert compute(term, net) == 2.0
@assert validate_term(term, net; verbose=false, rng=Xoshiro(1))
```

### When to Use This Pattern

- Geographic distance effects
- Pre-computed similarity scores
- Any dyad-level covariate stored as a matrix

### Adapting This Template

```julia
# Example: Binary covariate (same group membership)
struct SameGroup <: AbstractUserTerm
    membership::Vector{Int}
end

name(::SameGroup) = "samegroup"

function compute(t::SameGroup, net)
    count = 0.0
    for e in edges(net)
        i, j = src(e), dst(e)
        i <= length(t.membership) && j <= length(t.membership) || continue
        t.membership[i] == t.membership[j] && (count += 1.0)
    end
    return count
end

function change_stat(t::SameGroup, net, i::Int, j::Int)
    i <= length(t.membership) && j <= length(t.membership) || return 0.0
    return t.membership[i] == t.membership[j] ? 1.0 : 0.0
end
ERGM.is_dyad_dependent(::SameGroup) = false

@assert validate_term(SameGroup([1, 1, 2, 2]), net; verbose=false, rng=Xoshiro(1))
```

## InteractionTerm

Demonstrates combining multiple vertex attributes into an interaction effect.

### Implementation

```julia
struct InteractionTerm <: AbstractUserTerm
    attr1::Symbol
    attr2::Symbol
end

name(t::InteractionTerm) = "interact.$(t.attr1).$(t.attr2)"

# One vertex's value, 0.0 when it has none. `get_vertex_attribute(net, attr,
# v)` returns `nothing` when the attribute is absent or the vertex has no
# value; the zero-fill is reached only by raw compute calls, because the
# term DECLARES both attributes (bottom of this block) and an ERGMModel
# refuses a network on which either is absent or set on only some vertices
# — statnet refuses NA. The `::Float64` assertion keeps the arithmetic
# type-stable and allocation-free (the store is an untyped `Dict{Int,Any}`).
function _vertex_value(net, attr::Symbol, v::Int)
    x = get_vertex_attribute(net, attr, v)
    return x === nothing ? 0.0 : Float64(x)::Float64
end

# The per-edge contribution of (i,j): a_i b_j + a_j b_i
_interaction(t::InteractionTerm, net, i::Int, j::Int) =
    _vertex_value(net, t.attr1, i) * _vertex_value(net, t.attr2, j) +
    _vertex_value(net, t.attr1, j) * _vertex_value(net, t.attr2, i)

function compute(t::InteractionTerm, net)
    total = 0.0
    for e in edges(net)
        total += _interaction(t, net, Int(src(e)), Int(dst(e)))
    end
    return total
end

function change_stat(t::InteractionTerm, net, i::Int, j::Int)
    # Add-direction: the per-edge contribution of (i,j)
    return _interaction(t, net, i, j)
end

# The term is meaningless without its two attributes: declare them, and
# ERGMModel construction raises an ArgumentError naming a missing one
ERGM.required_vertex_attributes(t::InteractionTerm) = (t.attr1, t.attr2)
ERGM.is_dyad_dependent(::InteractionTerm) = false
```

### Usage

```julia
net = network(4; directed=true)
add_edge!(net, 1, 2)

# Set vertex attributes — on EVERY vertex: a declared attribute with a value
# on only some vertices is refused by ERGMModel (statnet refuses NA)
set_vertex_attribute!(net, :age, Dict(1 => 25.0, 2 => 30.0, 3 => 41.0, 4 => 52.0))
set_vertex_attribute!(net, :income, Dict(1 => 50000.0, 2 => 60000.0, 3 => 55000.0, 4 => 70000.0))

term = InteractionTerm(:age, :income)
compute(term, net)
# = age[1]*income[2] + age[2]*income[1]
# = 25*60000 + 30*50000 = 3000000.0
@assert compute(term, net) == 3000000.0
@assert validate_term(term, net; verbose=false, rng=Xoshiro(1))

# Misspell an attribute and the model, not the fit, tells you
err = try; ERGMModel(ERGMFormula([InteractionTerm(:age, :incom)]), net); nothing; catch e; e; end
@assert err isa ArgumentError && occursin(":incom", err.msg)
```

### When to Use This Pattern

- Cross-attribute effects (e.g., age × status)
- Sender-receiver attribute interactions
- Multi-dimensional homophily models

### Adapting This Template

```julia
# Example: Absolute difference interaction
struct DiffInteraction <: AbstractUserTerm
    attr1::Symbol
    attr2::Symbol
end

name(t::DiffInteraction) = "diffinteract.$(t.attr1).$(t.attr2)"

# Per-vertex read: 0 B and infers Float64 (the whole-Dict getter would
# allocate and infer Any on every sampler step — see the Benchmarking guide)
function _value(net, attr::Symbol, v::Int)
    x = get_vertex_attribute(net, attr, v)
    return x === nothing ? 0.0 : Float64(x)::Float64
end

# The per-edge contribution of (i, j): |Δattr1| · |Δattr2|
_diff_interaction(t::DiffInteraction, net, i::Int, j::Int) =
    abs(_value(net, t.attr1, i) - _value(net, t.attr1, j)) *
    abs(_value(net, t.attr2, i) - _value(net, t.attr2, j))

function compute(t::DiffInteraction, net)
    total = 0.0
    for e in edges(net)
        total += _diff_interaction(t, net, Int(src(e)), Int(dst(e)))
    end
    return total
end

change_stat(t::DiffInteraction, net, i::Int, j::Int) = _diff_interaction(t, net, i, j)
ERGM.required_vertex_attributes(t::DiffInteraction) = (t.attr1, t.attr2)
ERGM.is_dyad_dependent(::DiffInteraction) = false

@assert validate_term(DiffInteraction(:age, :income), net; verbose=false, rng=Xoshiro(1))
```

## Common Patterns Summary

### Counting Subgraphs

For terms that count specific subgraph patterns:

<!-- skip-check -->
```julia
function compute(::TriangleTerm, net)
    count = 0
    for i in vertices(net)
        for j in outneighbors(net, i)
            for k in outneighbors(net, j)
                k > i && has_edge(net, i, k) && (count += 1)
            end
        end
    end
    return count / 3  # Each triangle counted 3 times
end
```

### Using Attributes

For terms based on vertex properties, **declare** the attribute
(`ERGM.required_vertex_attributes(t::NodeMatchTerm) = (t.attr,)`). Inside an
`ERGMModel` a declared vertex attribute is then guaranteed to have a value on
every vertex — `ERGM._validate_formula` refuses a network on which it is
absent or partial, with `term '…' needs vertex attribute :attr on every
vertex … statnet refuses NA` — so the `get(attrs, v, nothing)` fallback
below is reached only by raw `compute` calls. An *undeclared* attribute read
is a [`validate_traits`](@ref) failure, because on a network lacking the
attribute the term would silently become an all-zero column (see
[Missing Attribute Values (NA)](@ref)):

<!-- skip-check -->
```julia
function compute(t::NodeMatchTerm, net)
    # get_vertex_attribute returns an empty Dict when the attribute is
    # absent; the `nothing` fallback only matters for raw compute calls —
    # a declared attribute is complete inside a model
    attrs = get_vertex_attribute(net, t.attr)
    count = 0
    for e in edges(net)
        get(attrs, src(e), nothing) == get(attrs, dst(e), nothing) && (count += 1)
    end
    return Float64(count)
end
ERGM.required_vertex_attributes(t::NodeMatchTerm) = (t.attr,)
```

### Conditional Statistics

For terms with conditional logic:

<!-- skip-check -->
```julia
function compute(t::ConditionalTerm, net)
    total = 0.0
    for e in edges(net)
        if meets_condition(net, src(e), dst(e), t.threshold)
            total += contribution(net, src(e), dst(e))
        end
    end
    return total
end
```

## Choosing a Template

| Your Term Needs | Start With |
|-----------------|------------|
| No parameters, simple counting | `ExampleTerm` |
| A tunable parameter | `TemplateTerm` |
| Edge-level data | `WeightedEdges` |
| An external covariate matrix | `DyadCovTerm` |
| Multiple vertex attributes | `InteractionTerm` |
