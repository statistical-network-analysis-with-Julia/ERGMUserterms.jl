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

The template implementations below re-define local copies of the bundled
terms so you can adapt them:

```julia
using ERGM, ERGMUserterms, Network
using Graphs: src, dst
import ERGM: name, compute, change_stat

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

# Edge attributes are keyed canonically: (i,j) directed, (min,max) undirected
_edge_key(net, i, j) = is_directed(net) ? (i, j) : minmax(i, j)

function compute(t::WeightedEdges, net)
    # get_edge_attribute returns an (empty) Dict even when the attribute is
    # absent, so sum over the network's edges with a default — this keeps
    # compute consistent with change_stat when the sampler adds unweighted
    # edges
    weights = get_edge_attribute(net, t.attr)
    total = 0.0
    for e in edges(net)
        total += Float64(get(weights, _edge_key(net, src(e), dst(e)), t.default))
    end
    return total
end

function change_stat(t::WeightedEdges, net, i::Int, j::Int)
    # Add-direction: the weight edge (i,j) carries (stored value, else the
    # default a fresh edge would get)
    weights = get_edge_attribute(net, t.attr)
    return Float64(get(weights, _edge_key(net, i, j), t.default))
end
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
    total = 0.0
    for e in edges(net)
        total += Float64(get(weights, _edge_key(net, src(e), dst(e)), 1.0))^2
    end
    return total
end

function change_stat(t::SquaredWeights, net, i::Int, j::Int)
    weights = get_edge_attribute(net, t.attr)
    w = Float64(get(weights, _edge_key(net, i, j), 1.0))
    return w^2
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
    # Undirected edges are stored canonically as (min, max); read the
    # covariate the same way compute() sees the edge
    if !is_directed(net)
        i, j = minmax(i, j)
    end
    i <= size(t.covariate, 1) && j <= size(t.covariate, 2) || return 0.0
    return t.covariate[i, j]
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

net = network(4; directed=true)
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
    return t.membership[i] == t.membership[j] ? 1.0 : 0.0
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
    # get_vertex_attribute returns an empty Dict (never nothing) when the
    # attribute is absent; missing values fall back to 0.0 via get()
    attrs1 = get_vertex_attribute(net, t.attr1)
    attrs2 = get_vertex_attribute(net, t.attr2)

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

    v1_i = get(attrs1, i, 0.0)
    v1_j = get(attrs1, j, 0.0)
    v2_i = get(attrs2, i, 0.0)
    v2_j = get(attrs2, j, 0.0)

    # Add-direction: the per-edge contribution of (i,j)
    return v1_i * v2_j + v1_j * v2_i
end
```

### Usage

```julia
net = network(4; directed=true)
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
    diff1 = abs(get(attrs1, i, 0.0) - get(attrs1, j, 0.0))
    diff2 = abs(get(attrs2, i, 0.0) - get(attrs2, j, 0.0))
    return diff1 * diff2
end
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

For terms based on vertex properties:

<!-- skip-check -->
```julia
function compute(t::NodeMatchTerm, net)
    # get_vertex_attribute returns an empty Dict when the attribute is
    # absent; vertices without a value fall back to `nothing` via get()
    attrs = get_vertex_attribute(net, t.attr)
    count = 0
    for e in edges(net)
        get(attrs, src(e), nothing) == get(attrs, dst(e), nothing) && (count += 1)
    end
    return Float64(count)
end
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
