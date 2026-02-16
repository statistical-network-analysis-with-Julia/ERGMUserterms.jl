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

## ExampleTerm

The simplest possible custom term. Counts edges weighted by vertex ID sum.

### Implementation

```julia
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
    has_edge(net, i, j) ? -(i + j) : Float64(i + j)
end
```

### Usage

```julia
term = ExampleTerm()
net = Network{Int}(; n=5, directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 3, 4)

compute(term, net)  # (1+2) + (3+4) = 10.0
change_stat(term, net, 2, 3)  # +(2+3) = 5.0 (edge absent)
change_stat(term, net, 1, 2)  # -(1+2) = -3.0 (edge present)
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
    i <= j && return 0.0  # This edge doesn't contribute
    return has_edge(net, i, j) ? -1.0 : 1.0
end
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
    # Change is just the parameter (positive if adding, negative if removing)
    return has_edge(net, i, j) ? -t.param : t.param
end
```

### Usage

```julia
# Create with different parameter types
term_float = TemplateTerm(2.5)
term_int = TemplateTerm(3)

net = Network{Int}(; n=5, directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)

compute(term_float, net)  # 2 * 2.5 = 5.0
change_stat(term_float, net, 3, 4)  # +2.5 (adding edge)
change_stat(term_float, net, 1, 2)  # -2.5 (removing edge)
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
    weight = t.decay ^ abs(i - j)
    return has_edge(net, i, j) ? -weight : weight
end
```

## WeightedEdges

Demonstrates accessing edge attributes from the network.

### Implementation

```julia
struct WeightedEdges <: AbstractUserTerm
    attr::Symbol
    WeightedEdges(attr::Symbol=:weight) = new(attr)
end

name(t::WeightedEdges) = "weightedges.$(t.attr)"

function compute(t::WeightedEdges, net)
    weights = get_edge_attribute(net, t.attr)
    isnothing(weights) && return Float64(ne(net))
    return Float64(sum(values(weights)))
end

function change_stat(t::WeightedEdges, net, i::Int, j::Int)
    weights = get_edge_attribute(net, t.attr)
    current_weight = if !isnothing(weights)
        get(weights, (i, j), 1.0)
    else
        1.0
    end
    return has_edge(net, i, j) ? -current_weight : current_weight
end
```

### Usage

```julia
net = Network{Int}(; n=5, directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)
set_edge_attribute!(net, :weight, (1,2), 3.0)
set_edge_attribute!(net, :weight, (2,3), 1.5)

term = WeightedEdges(:weight)
compute(term, net)  # 3.0 + 1.5 = 4.5
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

function compute(t::SquaredWeights, net)
    weights = get_edge_attribute(net, t.attr)
    isnothing(weights) && return Float64(ne(net))
    return sum(v^2 for v in values(weights))
end

function change_stat(t::SquaredWeights, net, i::Int, j::Int)
    weights = get_edge_attribute(net, t.attr)
    w = !isnothing(weights) ? get(weights, (i, j), 1.0) : 1.0
    return has_edge(net, i, j) ? -w^2 : w^2
end
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
    i <= size(t.covariate, 1) && j <= size(t.covariate, 2) || return 0.0
    val = t.covariate[i, j]
    return has_edge(net, i, j) ? -val : val
end
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

net = Network{Int}(; n=4, directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 3, 4)

term = DyadCovTerm(distance)
compute(term, net)  # distance[1,2] + distance[3,4] = 1.0 + 1.0 = 2.0
change_stat(term, net, 1, 3)  # +distance[1,3] = 2.0 (adding edge)
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
    t.membership[i] != t.membership[j] && return 0.0
    return has_edge(net, i, j) ? -1.0 : 1.0
end
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

function compute(t::InteractionTerm, net)
    attrs1 = get_vertex_attribute(net, t.attr1)
    attrs2 = get_vertex_attribute(net, t.attr2)
    (isnothing(attrs1) || isnothing(attrs2)) && return 0.0

    total = 0.0
    for e in edges(net)
        i, j = src(e), dst(e)
        v1_i = get(attrs1, i, 0.0)
        v1_j = get(attrs1, j, 0.0)
        v2_i = get(attrs2, i, 0.0)
        v2_j = get(attrs2, j, 0.0)
        total += v1_i * v2_j + v1_j * v2_i
    end
    return total
end

function change_stat(t::InteractionTerm, net, i::Int, j::Int)
    attrs1 = get_vertex_attribute(net, t.attr1)
    attrs2 = get_vertex_attribute(net, t.attr2)
    (isnothing(attrs1) || isnothing(attrs2)) && return 0.0

    v1_i = get(attrs1, i, 0.0)
    v1_j = get(attrs1, j, 0.0)
    v2_i = get(attrs2, i, 0.0)
    v2_j = get(attrs2, j, 0.0)

    delta = v1_i * v2_j + v1_j * v2_i
    return has_edge(net, i, j) ? -delta : delta
end
```

### Usage

```julia
net = Network{Int}(; n=4, directed=true)
add_edge!(net, 1, 2)

# Set vertex attributes
set_vertex_attribute!(net, :age, 1, 25.0)
set_vertex_attribute!(net, :age, 2, 30.0)
set_vertex_attribute!(net, :income, 1, 50000.0)
set_vertex_attribute!(net, :income, 2, 60000.0)

term = InteractionTerm(:age, :income)
compute(term, net)
# = age[1]*income[2] + age[2]*income[1]
# = 25*60000 + 30*50000 = 3000000.0
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

function compute(t::DiffInteraction, net)
    attrs1 = get_vertex_attribute(net, t.attr1)
    attrs2 = get_vertex_attribute(net, t.attr2)
    (isnothing(attrs1) || isnothing(attrs2)) && return 0.0
    total = 0.0
    for e in edges(net)
        i, j = src(e), dst(e)
        diff1 = abs(get(attrs1, i, 0.0) - get(attrs1, j, 0.0))
        diff2 = abs(get(attrs2, i, 0.0) - get(attrs2, j, 0.0))
        total += diff1 * diff2
    end
    return total
end

function change_stat(t::DiffInteraction, net, i::Int, j::Int)
    attrs1 = get_vertex_attribute(net, t.attr1)
    attrs2 = get_vertex_attribute(net, t.attr2)
    (isnothing(attrs1) || isnothing(attrs2)) && return 0.0
    diff1 = abs(get(attrs1, i, 0.0) - get(attrs1, j, 0.0))
    diff2 = abs(get(attrs2, i, 0.0) - get(attrs2, j, 0.0))
    delta = diff1 * diff2
    return has_edge(net, i, j) ? -delta : delta
end
```

## Common Patterns Summary

### Counting Subgraphs

For terms that count specific subgraph patterns:

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

For terms based on vertex properties:

```julia
function compute(t::NodeMatchTerm, net)
    attrs = get_vertex_attribute(net, t.attr)
    isnothing(attrs) && return 0.0
    count = 0
    for e in edges(net)
        get(attrs, src(e), nothing) == get(attrs, dst(e), nothing) && (count += 1)
    end
    return Float64(count)
end
```

### Conditional Statistics

For terms with conditional logic:

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
