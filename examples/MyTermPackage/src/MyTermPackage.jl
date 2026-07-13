"""
    MyTermPackage

A **package template** for shipping third-party ERGM terms.

Copy this directory, rename the module, and replace `ReciprocatedHomophily`
with your own term. It shows the complete contract a term must satisfy to be a
first-class citizen of ERGM.jl:

1. subtype `AbstractERGMTerm` (or `ERGMUserterms.AbstractUserTerm`);
2. implement the three interface methods — `name`, `compute`, `change_stat`
   (the last in the **add-direction** convention);
3. **declare its traits** through ERGM.jl's public term-trait protocol, so that
   `ERGMModel` construction validates the term exactly as it validates ERGM's
   own built-ins:

   | trait                            | this term          |
   |:---------------------------------|:-------------------|
   | `ERGM.required_vertex_attributes`| `(t.attr,)`        |
   | `ERGM.required_edge_attributes`  | `()` (the default) |
   | `ERGM.requires_directed`         | `true`             |
   | `ERGM.is_dyad_dependent`         | `true`             |
   | `Networks.supports_missing`      | `true`             |

`test/runtests.jl` asserts that `ERGMUserterms.validate_term` — which
exercises every one of those declarations — passes, and that the term fits
inside a real `ERGMModel`.
"""
module MyTermPackage

using ERGM
using ERGMUserterms
using Graphs
using Networks

# The interface generics are ERGM.jl's (they are really Networks.jl's, shared
# across the ecosystem). Import them by name so that your methods extend the
# one true function instead of defining a rival local one.
import ERGM: compute, change_stat, name

export ReciprocatedHomophily

"""
    ReciprocatedHomophily(attr::Symbol) <: AbstractUserTerm

Number of **mutual** dyads whose two vertices share the value of vertex
attribute `attr` — "reciprocated homophily". Unobserved (masked) dyads are
excluded: a pair with either arc masked contributes nothing, whatever its
stored face value.

# Example
```julia
net = network(5; directed=true)
set_vertex_attribute!(net, :group, Dict(v => isodd(v) ? "a" : "b" for v in 1:5))
add_edge!(net, 1, 3); add_edge!(net, 3, 1)     # mutual, both group "a"

compute(ReciprocatedHomophily(:group), net)    # 1.0
fit_ergm(net, [Edges(), ReciprocatedHomophily(:group)])
```
"""
struct ReciprocatedHomophily <: AbstractUserTerm
    attr::Symbol
end

name(t::ReciprocatedHomophily) = "recip_homophily.$(t.attr)"

# A dyad counts only if it is observed in both directions and its endpoints
# match on the attribute. Vertices without a value never match (matching the
# convention of ERGM's own nodematch).
function _counts(t::ReciprocatedHomophily, net, i::Int, j::Int)
    (is_missing_dyad(net, i, j) || is_missing_dyad(net, j, i)) && return false
    vals = get_vertex_attribute(net, t.attr)
    vi = get(vals, i, nothing)
    vj = get(vals, j, nothing)
    return !isnothing(vi) && isequal(vi, vj)
end

function compute(t::ReciprocatedHomophily, net)
    total = 0.0
    for e in edges(net)
        i, j = Int(src(e)), Int(dst(e))
        i < j || continue                       # count each mutual dyad once
        has_edge(net, j, i) || continue
        _counts(t, net, i, j) && (total += 1.0)
    end
    return total
end

# Add-direction change statistic: g(y⁺ᵢⱼ) − g(y⁻ᵢⱼ). Adding arc (i,j) completes
# a mutual dyad exactly when the reverse arc is already there. It must NOT look
# at whether (i,j) itself currently exists — that is the whole convention.
function change_stat(t::ReciprocatedHomophily, net, i::Int, j::Int)
    has_edge(net, j, i) || return 0.0
    return _counts(t, net, i, j) ? 1.0 : 0.0
end

# =============================================================================
# Trait declarations — the public protocol (ERGM.jl `src/terms/traits.jl`)
# =============================================================================

# 1. The vertex attribute the term reads. Without this declaration a model
#    built on a network lacking `:attr` would fit an all-zero column and
#    return a meaningless coefficient; with it, ERGMModel construction throws
#    the standard ArgumentError naming the attribute.
ERGM.required_vertex_attributes(t::ReciprocatedHomophily) = (t.attr,)

# 2. Mutuality is undefined on an undirected network, where every dyad is
#    trivially "reciprocated". Declaring the requirement makes ERGMModel
#    reject the term there instead of fitting a duplicate of nodematch.
ERGM.requires_directed(::ReciprocatedHomophily) = true

# 3. The change statistic reads another dyad (the reverse arc), so the term is
#    dyad-dependent: MPLE is not exact for it, and `show(::ERGMResult)` prints
#    the pseudo-likelihood caveat. (`true` is also ERGM's conservative default;
#    it is declared here to be explicit.)
ERGM.is_dyad_dependent(::ReciprocatedHomophily) = true

# 4. The statistic consults the missing-dyad mask (`_counts` above), so it is
#    invariant to the face value of any masked dyad. That — and only that — is
#    what `supports_missing` claims.
Networks.supports_missing(::ReciprocatedHomophily) = true

end # module
