"""
    ERGMUserterms.jl - Custom ERGM Term Development

Provides templates, utilities, and validation tools for developing custom ERGM terms.
Includes example terms and a comprehensive testing framework.

Port of the R ergm.userterms package from the StatNet collection.
"""
module ERGMUserterms

using ERGM
using Graphs
using Networks
using Random
using Statistics

# Extend the ERGM generics so user terms integrate with TermSet,
# summary_stats, and the estimation machinery
import ERGM: compute, change_stat, name

# Term development
export @ergm_term, validate_term, validate_traits, test_term
export AbstractUserTerm
# Re-export the interface functions users must extend
export compute, change_stat, name

# Templates and examples
export ExampleTerm, TemplateTerm
export WeightedEdges, DyadCovTerm, InteractionTerm

# Testing utilities
export change_stat_check, consistency_check
export benchmark_term, profile_term

# Documentation helpers
export term_signature, term_documentation

# =============================================================================
# Re-exported term interface (generics owned by ERGM.jl)
# =============================================================================
# The docstrings below document the interface from the term author's
# perspective; the generics themselves live in ERGM.jl.

"""
    name(term::AbstractERGMTerm) -> String

Return the term's descriptive, lowercase name, including any parameters
(e.g. `"template.2.0"`, `"weightedges.weight"`).

This generic is owned by ERGM.jl and re-exported here; add a method for
your term type. Without one, ERGM.jl falls back to `string(typeof(term))`.
"""
name

"""
    compute(term::AbstractERGMTerm, net) -> Float64

Compute the term's full network statistic `g(y)` on `net`.

This generic is owned by ERGM.jl and re-exported here; add a method for
your term type. It must be deterministic, must not modify the network, and
should handle empty and complete networks.
"""
compute

"""
    change_stat(term::AbstractERGMTerm, net, i::Int, j::Int) -> Float64

Compute the **add-direction** change statistic for dyad `(i, j)`:
`g(y⁺ᵢⱼ) − g(y⁻ᵢⱼ)`, i.e. the statistic with edge (i, j) present minus the
statistic with it absent, holding all other dyads fixed.

The value must **not** depend on whether the edge currently exists:
ERGM.jl's MPLE design matrix uses it directly, and the Metropolis–Hastings
sampler negates it for removal proposals. The toggle-direction idiom
`has_edge(net, i, j) ? -delta : delta` is wrong under this contract and is
rejected by [`change_stat_check`](@ref) / [`validate_term`](@ref).

This generic is owned by ERGM.jl and re-exported here; add a method for
your term type. For performance it should be O(degree), not O(edges).
"""
change_stat

# =============================================================================
# Term Development Infrastructure
# =============================================================================

"""
    AbstractUserTerm <: AbstractERGMTerm

Base type for user-defined ERGM terms.
Inherits from AbstractERGMTerm and provides additional validation hooks.
"""
abstract type AbstractUserTerm <: AbstractERGMTerm end

"""
    @ergm_term name body

Macro for defining custom ERGM terms. After evaluating `body`, it verifies
that the named type exists and that `compute`, `change_stat`, and `name`
methods are defined for it, warning about anything missing. (Numeric
consistency is checked separately with [`validate_term`](@ref) /
[`change_stat_check`](@ref), which need a term *instance* and a network.)

# Example
```julia
@ergm_term MyTerm begin
    struct MyTerm <: AbstractUserTerm
        param::Float64
    end

    name(t::MyTerm) = "myterm.\$(t.param)"

    function compute(t::MyTerm, net)
        # Implementation
    end

    function change_stat(t::MyTerm, net, i, j)
        # Implementation
    end
end
```
"""
macro ergm_term(termname, body)
    quote
        $(esc(body))
        _check_term_definition($(esc(termname)))
    end
end

# Does calling f on T hit only ERGM's AbstractERGMTerm fallback (i.e. the
# user supplied no method of their own)? hasmethod alone cannot tell,
# because ERGM defines error-throwing fallbacks for the whole interface.
function _only_fallback(f, argtypes::Type{<:Tuple})
    m = try
        which(f, argtypes)
    catch
        return true  # no method at all
    end
    return m.sig ==
           Tuple{typeof(f), AbstractERGMTerm, fieldtypes(argtypes)[2:end]...}
end

# Interface checks run by @ergm_term after the term is defined
function _check_term_definition(::Type{T}) where T
    if !(T <: AbstractERGMTerm)
        @warn "$(nameof(T)) does not subtype AbstractERGMTerm/AbstractUserTerm; " *
              "it will not work with ERGM.jl's TermSet"
        return nothing
    end
    _only_fallback(compute, Tuple{T, Any}) &&
        @warn "$(nameof(T)) has no compute(term, net) method"
    _only_fallback(change_stat, Tuple{T, Any, Int, Int}) &&
        @warn "$(nameof(T)) has no change_stat(term, net, i, j) method"
    _only_fallback(name, Tuple{T}) &&
        @warn "$(nameof(T)) has no name(term) method; ERGM will use the " *
              "default \"$(string(nameof(T)))\"-style name"
    return nothing
end
_check_term_definition(@nospecialize(x)) =
    @warn "@ergm_term body did not define a type (got $(typeof(x)))"

# =============================================================================
# Validation and Testing
# =============================================================================

"""
    validate_term(term::AbstractERGMTerm, net::Network; verbose=true, traits=true) -> Bool

Validate that a term is correctly implemented.

Checks the **method interface**:
- `compute()` returns a Real
- `change_stat()` returns correct values
- `name()` returns a non-empty string
- Change statistics are consistent with compute differences (add-direction
  convention, [`change_stat_check`](@ref))

and, unless `traits=false`, the term's **trait declarations** — see
[`validate_traits`](@ref) for exactly what is asserted about
`ERGM.required_vertex_attributes`, `ERGM.required_edge_attributes`,
`ERGM.requires_directed`/`requires_undirected`, `ERGM.is_dyad_dependent` and
`Networks.supports_missing`. The trait checks also assert that the term is
accepted by `ERGMModel` construction on `net`, so a term that passes
`validate_term` is a term ERGM.jl will fit.

`net` must carry every vertex/edge attribute the term declares (validate
against a network the term is meant for).
"""
function validate_term(term::AbstractERGMTerm, net::Network; verbose::Bool=true,
                       traits::Bool=true)
    valid = true

    # Check name
    try
        n = name(term)
        if isempty(n)
            verbose && @warn "Term name is empty"
            valid = false
        else
            verbose && @info "✓ name() returns: $n"
        end
    catch e
        verbose && @warn "name() failed: $e"
        valid = false
    end

    # Check compute
    try
        stat = compute(term, net)
        if !isa(stat, Real)
            verbose && @warn "compute() should return a Real, got $(typeof(stat))"
            valid = false
        else
            verbose && @info "✓ compute() returns: $stat ($(typeof(stat)))"
        end
    catch e
        verbose && @warn "compute() failed: $e"
        valid = false
    end

    # Validate change statistics
    n_vertices = nv(net)
    if n_vertices >= 2
        n_tests = min(10, n_vertices * (n_vertices - 1) ÷ 2)
        verbose && @info "Testing change_stat() with $n_tests random edges..."

        for _ in 1:n_tests
            i = rand(1:n_vertices)
            j = rand(1:(n_vertices - 1))
            j >= i && (j += 1)  # j ≠ i without discarding test draws

            try
                delta = change_stat(term, net, i, j)
                if !isa(delta, Real)
                    verbose && @warn "change_stat() should return a Real at ($i,$j)"
                    valid = false
                    break
                end
            catch e
                verbose && @warn "change_stat($i, $j) failed: $e"
                valid = false
                break
            end
        end

        if valid
            verbose && @info "✓ change_stat() returns valid values"
        end
    end

    # Consistency check
    if valid && n_vertices >= 2
        verbose && @info "Running consistency check..."
        consistent = change_stat_check(term, net; n_tests=5, verbose=false)
        if !consistent
            verbose && @warn "change_stat() values inconsistent with compute() differences"
            valid = false
        else
            verbose && @info "✓ change_stat() consistent with compute()"
        end
    end

    # Trait declarations (the public protocol in ERGM.jl `src/terms/traits.jl`)
    if valid && traits
        valid &= validate_traits(term, net; verbose=verbose)
    end

    return valid
end

# =============================================================================
# Trait validation
# =============================================================================
#
# A term does not only compute a number: it *declares* what it needs (vertex
# and edge attributes), where it is defined (directed / undirected networks),
# what it depends on (other dyads), and what it does with unobserved dyads.
# ERGM.jl acts on those declarations at model construction, so a wrong
# declaration is as damaging as a wrong change statistic — and, unlike a wrong
# change statistic, nothing else in the stack will catch it. These checks are
# the harness for the declarations.

# The term's observable behaviour on a network: the full statistic plus every
# add-direction change statistic. Two networks agreeing here are
# indistinguishable to ERGM.jl's estimation and sampling machinery.
function _term_fingerprint(term::AbstractERGMTerm, net::Network)
    n = Int(nv(net))
    vals = Float64[Float64(compute(term, net))]
    for i in 1:n, j in 1:n
        i == j && continue
        (!is_directed(net) && j < i) && continue
        push!(vals, Float64(change_stat(term, net, i, j)))
    end
    return vals
end

_fingerprints_differ(a, b; tol=1e-10) =
    length(a) != length(b) || any(abs.(a .- b) .> tol)

# Two type-preserving perturbations of an attribute's values. Each is applied
# to a copy; a term "reads" the attribute if either changes its fingerprint.
#
# - `:rotate` maps every value to the next distinct value (detects terms whose
#   statistic depends on the values themselves, e.g. nodecov, absdiff)
# - `:collapse` maps every value to a single one (detects terms that depend
#   only on the *partition* the values induce, e.g. nodematch — which `:rotate`
#   leaves invariant)
#
# Returns `false` when the attribute cannot be perturbed (absent, or constant
# — then it carries no information and dependence is undetectable).
function _perturb!(values::Dict, setter!, mode::Symbol)
    isempty(values) && return false
    uniq = unique(collect(Base.values(values)))
    length(uniq) < 2 && return false
    if mode === :rotate
        nxt = Dict(uniq[k] => uniq[mod1(k + 1, length(uniq))] for k in eachindex(uniq))
        for (key, val) in collect(values)
            setter!(key, nxt[val])
        end
    else  # :collapse
        for (key, _) in collect(values)
            setter!(key, first(uniq))
        end
    end
    return true
end

# Does `term`'s behaviour depend on vertex/edge attribute `attr`?
# `missing` when the attribute cannot be perturbed (dependence undecidable).
function _reads_attribute(term::AbstractERGMTerm, net::Network, attr::Symbol,
                          kind::Symbol)
    base = _term_fingerprint(term, net)
    decided = false
    for mode in (:rotate, :collapse)
        work = deepcopy(net)
        vals = kind === :vertex ? get_vertex_attribute(work, attr) :
                                  get_edge_attribute(work, attr)
        setter! = kind === :vertex ?
            (v, val) -> set_vertex_attribute!(work, attr, v, val) :
            ((i, j), val) -> set_edge_attribute!(work, attr, i, j, val)
        _perturb!(vals, setter!, mode) || continue
        decided = true
        # A term that errors on perturbed values is certainly reading them
        moved = try
            _fingerprints_differ(base, _term_fingerprint(term, work))
        catch
            true
        end
        moved && return true
    end
    return decided ? false : missing
end

# An undirected/directed structural twin of `net`, carrying the same vertex and
# edge attributes, used to check that a direction requirement is enforced.
function _direction_twin(net::Network, directed::Bool)
    twin = network(Int(nv(net)); directed=directed)
    for e in edges(net)
        add_edge!(twin, src(e), dst(e))
    end
    for attr in list_vertex_attributes(net)
        for (v, val) in get_vertex_attribute(net, attr)
            set_vertex_attribute!(twin, attr, v, val)
        end
    end
    for attr in list_edge_attributes(net)
        for ((i, j), val) in get_edge_attribute(net, attr)
            has_edge(twin, i, j) && set_edge_attribute!(twin, attr, i, j, val)
        end
    end
    return twin
end

_throws_argument_error(f) = try
    f()
    false
catch e
    e isa ArgumentError
end

"""
    validate_traits(term::AbstractERGMTerm, net::Network; verbose=true) -> Bool

Exercise a term's declarations under ERGM.jl's public term-trait protocol
(`ERGM.required_vertex_attributes`, `ERGM.required_edge_attributes`,
`ERGM.requires_directed`, `ERGM.requires_undirected`,
`ERGM.is_dyad_dependent`, `Networks.supports_missing`) against `net`.

Checked, in order:

1. **Attributes are the ones it reads.** Every attribute the term declares
   must exist on `net` and must actually move the term's statistic when its
   values are perturbed (otherwise the declaration is a warning: an
   over-declaration needlessly narrows the networks the term can be fitted
   to). Conversely, a term whose statistic moves when an **undeclared vertex
   attribute** is perturbed fails: it reads an attribute ERGM.jl will not
   validate, so on a network lacking that attribute it would silently become
   an all-zero design column instead of an error. (An undeclared *edge*
   attribute is only warned about — terms like [`WeightedEdges`](@ref) read
   one but fall back to a default, so they genuinely do not require it.)
2. **Direction requirements are honoured.** A term declaring
   `requires_directed` must be rejected by `ERGMModel` construction on the
   undirected twin of `net` (and symmetrically for `requires_undirected`).
3. **Dyad-independence claims hold.** A term declaring
   `is_dyad_dependent(term) == false` must return the same change statistic
   for a dyad no matter how *other* dyads are toggled. (`true`, the
   conservative default, asserts nothing.)
4. **Missing-data claims hold.** A term declaring
   `Networks.supports_missing(term) == true` must return the same statistic
   however the face value of a **masked** dyad is flipped — that is what it
   means to honour the mask. (`false`, the default, asserts nothing: the
   statistic then counts masked dyads at face value, and the missing-data
   treatment lives in the estimator.)

Finally, `ERGMModel` construction with the term on `net` must succeed.
"""
function validate_traits(term::AbstractERGMTerm, net::Network; verbose::Bool=true)
    valid = true
    n = Int(nv(net))

    # --- 1. Attribute declarations ------------------------------------------
    for (kind, declared, present) in
            ((:vertex, ERGM.required_vertex_attributes(term), list_vertex_attributes(net)),
             (:edge, ERGM.required_edge_attributes(term), list_edge_attributes(net)))
        for attr in declared
            if !(attr in present)
                verbose && @warn "term declares required $kind attribute :$attr, " *
                                 "which the validation network does not have; " *
                                 "ERGMModel construction would reject it"
                valid = false
            end
        end
        for attr in present
            reads = _reads_attribute(term, net, attr, kind)
            reads === missing && continue      # constant attribute: undecidable
            if reads && !(attr in declared)
                msg = "term's statistic depends on $kind attribute :$attr but " *
                      "does not declare it (ERGM.required_$(kind)_attributes)"
                if kind === :vertex
                    verbose && @warn msg * "; on a network without :$attr it would " *
                                     "silently evaluate to zero instead of raising"
                    valid = false
                else
                    verbose && @warn msg * "; declare it unless the term has a " *
                                     "meaningful default when the attribute is absent"
                end
            elseif !reads && attr in declared
                verbose && @warn "term declares required $kind attribute :$attr but " *
                                 "its statistic does not depend on it on this network " *
                                 "(over-declaration needlessly restricts the term)"
            end
        end
    end
    if valid && verbose
        @info "✓ declared attributes match the ones the term reads"
    end

    # --- 2. Direction requirements ------------------------------------------
    req_dir = ERGM.requires_directed(term)
    req_undir = ERGM.requires_undirected(term)
    if req_dir && req_undir
        verbose && @warn "term declares both requires_directed and " *
                         "requires_undirected: no network can satisfy it"
        valid = false
    elseif req_dir && !is_directed(net)
        verbose && @warn "term declares requires_directed but was validated on an " *
                         "undirected network; validate it on a directed one"
        valid = false
    elseif req_undir && is_directed(net)
        verbose && @warn "term declares requires_undirected but was validated on a " *
                         "directed network; validate it on an undirected one"
        valid = false
    elseif req_dir || req_undir
        twin = _direction_twin(net, !is_directed(net))
        if _throws_argument_error(() -> ERGMModel(ERGMFormula([term]), twin))
            verbose && @info "✓ direction requirement enforced by ERGMModel " *
                             "on a $(is_directed(twin) ? "directed" : "undirected") network"
        else
            verbose && @warn "term declares requires_$(req_dir ? "directed" : "undirected") " *
                             "but ERGMModel construction accepted it on a " *
                             "$(is_directed(twin) ? "directed" : "undirected") network"
            valid = false
        end
    end

    # --- 3. Dyad-dependence declaration -------------------------------------
    # Every *other* dyad is a candidate: a dyad-dependent change statistic may
    # read the reverse arc (mutual), an incident dyad (triangle, star), or a
    # distant one (geodesic terms), so sampling toggles at random would miss
    # most violations. All other dyads are toggled for each of a few source
    # dyads (capped on big networks).
    if !ERGM.is_dyad_dependent(term) && n >= 3
        work = deepcopy(net)
        sources = unique(_random_dyad(n) for _ in 1:4)
        for (i, j) in sources
            baseline = change_stat(term, work, i, j)
            for (k, l) in _other_dyads(work, i, j)
                saved = _edge_attr_snapshot(work, k, l)
                had = has_edge(work, k, l)
                had ? rem_edge!(work, k, l) : add_edge!(work, k, l)
                _edge_attr_restore!(work, k, l, saved)
                toggled = change_stat(term, work, i, j)
                had ? add_edge!(work, k, l) : rem_edge!(work, k, l)
                _edge_attr_restore!(work, k, l, saved)

                if abs(baseline - toggled) > 1e-10
                    verbose && @warn "term declares is_dyad_dependent = false, but " *
                                     "its change statistic at dyad ($i,$j) moved from " *
                                     "$baseline to $toggled when dyad ($k,$l) was " *
                                     "toggled: it IS dyad-dependent"
                    valid = false
                    break
                end
            end
            valid || break
        end
        valid && verbose && @info "✓ is_dyad_dependent = false holds under toggling " *
                                  "of other dyads"
    end

    # --- 4. Missing-data declaration ----------------------------------------
    if supports_missing(term) && n >= 2
        work = deepcopy(net)
        for _ in 1:5
            i, j = _random_dyad(n)
            is_missing_dyad(work, i, j) && continue
            set_missing_dyad!(work, i, j)
            saved = _edge_attr_snapshot(work, i, j)
            had = has_edge(work, i, j)
            s0 = compute(term, work)
            had ? rem_edge!(work, i, j) : add_edge!(work, i, j)
            _edge_attr_restore!(work, i, j, saved)
            s1 = compute(term, work)
            had ? add_edge!(work, i, j) : rem_edge!(work, i, j)
            _edge_attr_restore!(work, i, j, saved)
            delete_missing_dyad!(work, i, j)

            if abs(s0 - s1) > 1e-10
                verbose && @warn "term declares supports_missing = true, but its " *
                                 "statistic changed ($s0 → $s1) when the face value " *
                                 "of the MASKED dyad ($i,$j) was flipped: it counts " *
                                 "unobserved dyads at face value"
                valid = false
                break
            end
        end
        valid && verbose && @info "✓ supports_missing = true holds: masked dyads do " *
                                  "not enter the statistic"
    end

    # --- 5. ERGM model construction accepts the term -------------------------
    try
        ERGMModel(ERGMFormula([Edges(), term]), net)
        verbose && @info "✓ accepted by ERGMModel construction"
    catch e
        verbose && @warn "ERGMModel construction rejected the term on this network: $e"
        valid = false
    end

    return valid
end

# Every dyad of `net` other than `(i, j)` itself (and, on an undirected
# network, its mirror). Capped at 250 randomly chosen dyads so the check stays
# cheap on large networks.
function _other_dyads(net::Network, i::Int, j::Int; cap::Int=250)
    n = Int(nv(net))
    directed = is_directed(net)
    out = Tuple{Int, Int}[]
    for k in 1:n, l in 1:n
        k == l && continue
        (!directed && l < k) && continue
        (k, l) == (i, j) && continue
        (!directed && (k, l) == minmax(i, j)) && continue
        push!(out, (k, l))
    end
    length(out) > cap && (out = out[randperm(length(out))[1:cap]])
    return out
end

# Uniform draw from the off-diagonal dyads, without discarding draws
function _random_dyad(n::Int)
    i = rand(1:n)
    j = rand(1:(n - 1))
    j >= i && (j += 1)
    return i, j
end

# Edge attributes are treated as exogenous dyadic data by the harness:
# Networks.rem_edge! deletes them, so they are snapshotted before a toggle
# and restored afterwards. Without this, any term reading edge attributes
# would spuriously fail the brute-force comparison.
function _edge_attr_snapshot(net::Network, i::Int, j::Int)
    key = is_directed(net) ? (i, j) : minmax(i, j)
    return [attr => attrs[key] for (attr, attrs) in net.edge_attrs
            if haskey(attrs, key)]
end

function _edge_attr_restore!(net::Network, i::Int, j::Int, saved)
    key = is_directed(net) ? (i, j) : minmax(i, j)
    for (attr, val) in saved
        net.edge_attrs[attr][key] = val
    end
    return net
end

# Brute-force add-direction change statistic: g(y⁺ij) − g(y⁻ij) computed
# by actually toggling the dyad. Restores the network (including the
# dyad's edge attributes) to its original state.
function _brute_change_stat(term::AbstractERGMTerm, net::Network, i::Int, j::Int)
    had = has_edge(net, i, j)
    saved = _edge_attr_snapshot(net, i, j)

    had && rem_edge!(net, i, j)
    s0 = compute(term, net)
    add_edge!(net, i, j)
    _edge_attr_restore!(net, i, j, saved)
    s1 = compute(term, net)
    if !had
        rem_edge!(net, i, j)
        _edge_attr_restore!(net, i, j, saved)
    end
    return s1 - s0
end

# Check the add-direction invariant at dyad (i,j): change_stat must equal
# compute(edge present) − compute(edge absent) and must not depend on the
# dyad's current state
function _dyad_consistent(term::AbstractERGMTerm, test_net::Network,
                          i::Int, j::Int, tol::Float64, verbose::Bool)
    predicted = change_stat(term, test_net, i, j)
    expected = _brute_change_stat(term, test_net, i, j)

    if abs(predicted - expected) > tol
        verbose && @warn "Inconsistency at dyad ($i,$j): " *
                         "change_stat=$predicted, brute-force=$expected"
        return false
    end

    # State independence: same value with the dyad toggled
    had = has_edge(test_net, i, j)
    saved = _edge_attr_snapshot(test_net, i, j)
    had ? rem_edge!(test_net, i, j) : add_edge!(test_net, i, j)
    _edge_attr_restore!(test_net, i, j, saved)
    predicted_toggled = change_stat(term, test_net, i, j)
    had ? add_edge!(test_net, i, j) : rem_edge!(test_net, i, j)
    _edge_attr_restore!(test_net, i, j, saved)

    if abs(predicted_toggled - expected) > tol
        verbose && @warn "State-dependent change_stat at dyad ($i,$j): " *
                         "$predicted with edge state as-is vs " *
                         "$predicted_toggled after toggling. change_stat " *
                         "must return the add-direction change regardless " *
                         "of whether the edge exists"
        return false
    end

    return true
end

"""
    change_stat_check(term::AbstractERGMTerm, net::Network; n_tests=10, verbose=true) -> Bool

Verify that `change_stat` returns the **add-direction** change statistic
required by ERGM.jl: `compute` with edge (i,j) present minus `compute` with
it absent, independent of the dyad's current state. Random dyads are
toggled and restored.
"""
function change_stat_check(term::AbstractERGMTerm, net::Network;
                           n_tests::Int=10, verbose::Bool=true, tol::Float64=1e-10)
    test_net = deepcopy(net)
    n = nv(test_net)
    n < 2 && return true
    all_passed = true

    for _ in 1:n_tests
        i = rand(1:n)
        j = rand(1:(n - 1))
        j >= i && (j += 1)  # uniform over j ≠ i, without discarding draws

        all_passed &= _dyad_consistent(term, test_net, i, j, tol, verbose)
    end

    return all_passed
end

"""
    consistency_check(term::AbstractERGMTerm, net::Network; exhaustive=false) -> Bool

Full consistency check between `compute` and the add-direction
`change_stat` (see [`change_stat_check`](@ref)).
If exhaustive=true, checks all possible dyads (can be slow for large networks).
"""
function consistency_check(term::AbstractERGMTerm, net::Network;
                           exhaustive::Bool=false, tol::Float64=1e-10)
    n = nv(net)
    n < 2 && return true
    test_net = deepcopy(net)

    dyads_to_check = if exhaustive
        [(i, j) for i in 1:n for j in 1:n if i != j]
    else
        unique((rand(1:n), rand(1:n)) for _ in 1:min(100, n * n))
    end

    for (i, j) in dyads_to_check
        i == j && continue
        _dyad_consistent(term, test_net, i, j, tol, false) || return false
    end

    return true
end

"""
    test_term(term::AbstractERGMTerm; n_vertices=20, density=0.1, n_tests=100) -> Bool

Comprehensive test suite for an ERGM term.
"""
function test_term(term::AbstractERGMTerm;
                   n_vertices::Int=20,
                   density::Float64=0.1,
                   n_tests::Int=100)

    # Create test network
    net = Network{Int}(; n=n_vertices, directed=true)
    n_edges = round(Int, density * n_vertices * (n_vertices - 1))
    for _ in 1:n_edges
        i, j = rand(1:n_vertices), rand(1:n_vertices)
        i != j && add_edge!(net, i, j)
    end

    println("Testing term: $(name(term))")
    println("=" ^ 50)

    # Run validation
    valid = validate_term(term, net; verbose=true)

    # Additional tests
    if valid
        println("\nAdditional tests:")

        # Test on empty network
        empty_net = Network{Int}(; n=n_vertices)
        try
            stat_empty = compute(term, empty_net)
            println("✓ Works on empty network: $stat_empty")
        catch e
            println("✗ Failed on empty network: $e")
            valid = false
        end

        # Test on complete network
        complete_net = Network{Int}(; n=min(10, n_vertices))
        for i in 1:nv(complete_net), j in 1:nv(complete_net)
            i != j && add_edge!(complete_net, i, j)
        end
        try
            stat_complete = compute(term, complete_net)
            println("✓ Works on complete network: $stat_complete")
        catch e
            println("✗ Failed on complete network: $e")
            valid = false
        end
    end

    println("=" ^ 50)
    println(valid ? "All tests PASSED" : "Some tests FAILED")

    return valid
end

"""
    benchmark_term(term::AbstractERGMTerm, net::Network; n_iter=1000) -> NamedTuple

Benchmark compute() and change_stat() performance.
"""
function benchmark_term(term::AbstractERGMTerm, net::Network; n_iter::Int=1000)
    # Benchmark compute
    compute_times = Float64[]
    for _ in 1:n_iter
        t = @elapsed compute(term, net)
        push!(compute_times, t)
    end

    # Benchmark change_stat
    change_times = Float64[]
    n = nv(net)
    for _ in 1:n_iter
        i = rand(1:n)
        j = rand(1:(n - 1))
        j >= i && (j += 1)
        t = @elapsed change_stat(term, net, i, j)
        push!(change_times, t)
    end

    return (
        compute_mean = mean(compute_times),
        compute_std = std(compute_times),
        change_stat_mean = mean(change_times),
        change_stat_std = std(change_times),
        speedup = mean(compute_times) / mean(change_times)
    )
end

"""
    profile_term(term::AbstractERGMTerm; sizes=[10, 20, 40], density=0.1,
                 n_iter=200) -> Vector{NamedTuple}

Profile a term's `compute`/`change_stat` cost across network sizes.
Returns one `benchmark_term` result per size (with the size added), useful
for spotting change statistics that do not scale (e.g. O(edges) instead of
O(degree)).
"""
function profile_term(term::AbstractERGMTerm;
                      sizes::Vector{Int}=[10, 20, 40],
                      density::Float64=0.1,
                      n_iter::Int=200)
    results = NamedTuple[]
    for n in sizes
        net = Network{Int}(; n=n, directed=true)
        n_edges = round(Int, density * n * (n - 1))
        for _ in 1:n_edges
            i, j = rand(1:n), rand(1:n)
            i != j && add_edge!(net, i, j)
        end
        b = benchmark_term(term, net; n_iter=n_iter)
        push!(results, (n_vertices=n, ne=ne(net), b...))
    end
    return results
end

# =============================================================================
# Example Terms and Templates
# =============================================================================

"""
    ExampleTerm <: AbstractUserTerm

Example custom term that counts edges weighted by vertex ID sum.
Use this as a template for creating your own terms.

All `change_stat` methods must return the **add-direction** change
statistic `g(y⁺ij) − g(y⁻ij)` — the statistic with edge (i,j) present minus
with it absent — independent of whether the edge currently exists. ERGM.jl's
MPLE design matrix and MH sampler both require this convention.
"""
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
    # Adding edge (i,j) adds i + j to the statistic
    return Float64(i + j)
end

"""
    TemplateTerm{T} <: AbstractUserTerm

A parameterized template term. Copy and modify this for your own terms.

# Fields
- `param::T`: A parameter value
- `attr::Symbol`: An attribute name (optional)
"""
struct TemplateTerm{T} <: AbstractUserTerm
    param::T
    attr::Symbol

    TemplateTerm(param::T; attr::Symbol=:none) where T = new{T}(param, attr)
end

name(t::TemplateTerm) = "template.$(t.param)"

function compute(t::TemplateTerm, net)
    # Template: count edges multiplied by parameter
    return Float64(ne(net)) * t.param
end

function change_stat(t::TemplateTerm, net, i::Int, j::Int)
    # Template: adding edge (i,j) increases the statistic by the parameter
    # (add-direction, state-independent)
    return Float64(t.param)
end

"""
    WeightedEdges <: AbstractUserTerm

Sum of edge weights, with `default` (1.0) for edges lacking the weight
attribute. Demonstrates accessing edge attributes.

Note `get_edge_attribute` returns an (empty) Dict even when the attribute
is absent, so `compute` sums over the network's edges with a default rather
than over the attribute dictionary — that keeps it consistent with
`change_stat` when edges are added without weights (e.g. by the MH sampler).
"""
struct WeightedEdges <: AbstractUserTerm
    attr::Symbol
    default::Float64
    WeightedEdges(attr::Symbol=:weight; default::Float64=1.0) = new(attr, default)
end

name(t::WeightedEdges) = "weightedges.$(t.attr)"

# Edge attributes are keyed canonically: (i,j) for directed networks,
# (min,max) for undirected
_edge_key(net, i::Int, j::Int) = is_directed(net) ? (i, j) : minmax(i, j)

function compute(t::WeightedEdges, net)
    weights = get_edge_attribute(net, t.attr)
    total = 0.0
    for e in edges(net)
        total += Float64(get(weights, _edge_key(net, src(e), dst(e)), t.default))
    end
    return total
end

function change_stat(t::WeightedEdges, net, i::Int, j::Int)
    # Add-direction: the weight edge (i,j) carries — its stored attribute
    # if present, otherwise the default a fresh edge would get
    weights = get_edge_attribute(net, t.attr)
    return Float64(get(weights, _edge_key(net, i, j), t.default))
end

"""
    DyadCovTerm <: AbstractUserTerm

Dyadic covariate term: sum over edges of covariate values.
Demonstrates using a matrix covariate.
"""
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
    # covariate the same way compute() will see the edge
    if !is_directed(net)
        i, j = minmax(i, j)
    end
    i <= size(t.covariate, 1) && j <= size(t.covariate, 2) || return 0.0
    return t.covariate[i, j]
end

"""
    InteractionTerm <: AbstractUserTerm

Interaction between two node attributes.
"""
struct InteractionTerm <: AbstractUserTerm
    attr1::Symbol
    attr2::Symbol
end

name(t::InteractionTerm) = "interact.$(t.attr1).$(t.attr2)"

function compute(t::InteractionTerm, net)
    # Note: get_vertex_attribute returns an empty Dict (never nothing) when
    # the attribute is absent; missing values fall back to 0.0 via get()
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

    # Add-direction change: the per-edge contribution of (i,j)
    return v1_i * v2_j + v1_j * v2_i
end

# =============================================================================
# Trait declarations for the bundled terms
# =============================================================================
#
# The public term-trait protocol (ERGM.jl `src/terms/traits.jl`) is how a term
# tells ERGM.jl what it needs and where it belongs. Every term you write should
# declare the traits that are not at their default — see the worked example in
# `examples/MyTermPackage/`.

# All bundled example terms are covariate-only: their change statistics
# never read the state of other dyads. Declaring dyad-independence opts
# them out of ERGM.jl's conservative fallback (`is_dyad_dependent = true`
# for unknown term types), which otherwise triggers the pseudo-likelihood
# caveat in `show` and a conservative MCMLE bridge reference. User terms
# should do the same when (and only when) they are covariate-only.
ERGM.is_dyad_dependent(::ExampleTerm) = false
ERGM.is_dyad_dependent(::TemplateTerm) = false
ERGM.is_dyad_dependent(::WeightedEdges) = false
ERGM.is_dyad_dependent(::DyadCovTerm) = false
ERGM.is_dyad_dependent(::InteractionTerm) = false

# `InteractionTerm` reads two vertex attributes and has no meaning without
# them: on a network missing either, `get_vertex_attribute` returns an empty
# Dict and the statistic collapses to a constant zero — an all-zero design
# column and a meaningless coefficient. Declaring them makes ERGMModel
# construction raise the standard ArgumentError instead.
ERGM.required_vertex_attributes(t::InteractionTerm) = (t.attr1, t.attr2)

# `WeightedEdges` deliberately declares NO required edge attribute: an edge
# without the weight attribute (e.g. one the MH sampler just proposed) counts
# at `t.default`, so the term is perfectly well defined on a network that has
# never heard of `t.attr`. Declaring it would reject those networks. Required
# means "an error if absent", not "read if present".

# =============================================================================
# Documentation Helpers
# =============================================================================

"""
    term_signature(term::AbstractERGMTerm) -> String

Generate a signature string for a term showing its fields and types.
"""
function term_signature(term::AbstractERGMTerm)
    T = typeof(term)
    fields = fieldnames(T)
    types = [fieldtype(T, f) for f in fields]

    sig = "$(nameof(T))"
    if !isempty(fields)
        params = ["$f::$t" for (f, t) in zip(fields, types)]
        sig *= "(" * join(params, ", ") * ")"
    end
    return sig
end

"""
    term_documentation(term::AbstractERGMTerm) -> String

Generate documentation for a term.
"""
function term_documentation(term::AbstractERGMTerm)
    T = typeof(term)
    doc = """
    # $(nameof(T))

    Signature: $(term_signature(term))
    Name: $(name(term))

    ## Fields
    """

    for f in fieldnames(T)
        doc *= "- `$f::$(fieldtype(T, f))` = $(getfield(term, f))\n"
    end

    doc *= """

    ## Interface
    - `name(term)` -> String
    - `compute(term, net)` -> Float64
    - `change_stat(term, net, i, j)` -> Float64
    """

    return doc
end

end # module
