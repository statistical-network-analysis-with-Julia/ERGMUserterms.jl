"""
    ERGMUserterms.jl - Custom ERGM Term Development

Provides templates, utilities, and validation tools for developing custom ERGM terms.
Includes example terms and a comprehensive testing framework.

Port of the R ergm.userterms package from the StatNet collection. What is
ported is the *role* of that package — a validated, copyable starting point
for third-party terms — not its contents: R's C-level `changestats.users`
template and its `mindegree` example term are not ported (the Julia
counterpart of the C template is `examples/MyTermPackage/`; a Julia term is
plain Julia, so there is no C skeleton to fill in).

The interface generics a term extends — `name`, `compute`, `change_stat` —
are ERGM.jl's (`compute`/`name` are really Networks.jl's shared statistic
protocol); this package imports and re-exports them, so `import ERGM: …` and
`import ERGMUserterms: …` name the same functions. Edge endpoints `src`/`dst`
are `Graphs.src`/`Graphs.dst`, which neither Networks.jl nor ERGM.jl
re-exports: a term iterating `edges(net)` needs `using Graphs: src, dst`.
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

This generic is owned by ERGM.jl (it is `Networks.name`, the ecosystem's
shared statistic protocol) and re-exported here: `ERGMUserterms.name ===
ERGM.name === Networks.name`. Add a method for your term type after
importing it by name — `import ERGM: name` or `import ERGMUserterms: name`,
the same function; a bare `using` followed by `name(::MyTerm) = …` defines a
*local* `name` that ERGM.jl never calls. Without a method, ERGM.jl falls
back to `string(typeof(term))`.

Use R ergm's label where one exists (`"gwdeg.fixed.0.5"`, `"nodematch.group"`)
so coefficient tables line up with statnet output. ERGM.jl also has a
two-argument `name(term, net)` for labels that depend on the network's
direction (`"gwesp.OTP.fixed.0.5"` on a digraph); it falls back to
`name(term)`, so a user term needs only the one-argument method.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
import ERGM: name, compute, change_stat     # extend the shared generics
struct OneEdge <: AbstractUserTerm end
name(::OneEdge) = "oneedge"
compute(::OneEdge, net) = Float64(ne(net))
change_stat(::OneEdge, net, i::Int, j::Int) = 1.0
name(OneEdge())                                   # "oneedge"
name(OneEdge()) == ERGM.name(OneEdge())           # true — one generic, not a local copy
validate_term(OneEdge(), network(4; directed=true); verbose=false, rng=Xoshiro(1))   # true
```
"""
name

"""
    compute(term::AbstractERGMTerm, net) -> Float64

Compute the term's full network statistic `g(y)` on `net`.

This generic is owned by ERGM.jl (it is `Networks.compute`, the ecosystem's
shared statistic protocol) and re-exported here: `ERGMUserterms.compute ===
ERGM.compute === Networks.compute`. Add a method for your term type after
`import ERGM: compute` (or `import ERGMUserterms: compute` — the same
function). It must be deterministic, must not modify the network, and should
handle empty and complete networks.

Edge endpoints are `Graphs.src`/`Graphs.dst`, which Networks.jl does not
re-export — a `compute` that iterates `edges(net)` needs `using Graphs: src,
dst`. A vertex attribute the term *declares* (`ERGM.required_vertex_attributes`)
is guaranteed complete inside an `ERGMModel` (ERGM refuses a partial one, as
statnet refuses NA), so `get(attrs, v, default)` fallbacks only ever fire on
raw `compute` calls.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
using Graphs: src, dst
import ERGM: name, compute, change_stat
struct Descending <: AbstractUserTerm end      # arcs i -> j with i > j
name(::Descending) = "descending"
compute(::Descending, net) = Float64(count(src(e) > dst(e) for e in edges(net)))
change_stat(::Descending, net, i::Int, j::Int) = i > j ? 1.0 : 0.0
net = network(4; directed=true)
add_edge!(net, 3, 1); add_edge!(net, 1, 2)
compute(Descending(), net)                                            # 1.0
validate_term(Descending(), net; verbose=false, rng=Xoshiro(1))       # true
```
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

This generic is owned by ERGM.jl and re-exported here
(`ERGMUserterms.change_stat === ERGM.change_stat`); add a method for your
term type after `import ERGM: change_stat`. For performance it should be
O(degree), not O(edges).

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
import ERGM: name, compute, change_stat
struct Mutuals <: AbstractUserTerm end          # mutual dyads, counted once
name(::Mutuals) = "mutuals"
compute(::Mutuals, net) =
    Float64(count(has_edge(net, i, j) && has_edge(net, j, i)
                  for i in 1:nv(net) for j in 1:nv(net) if i < j))
# Adding i -> j completes a mutual dyad iff j -> i is there; the value does
# not consult has_edge(net, i, j) itself
change_stat(::Mutuals, net, i::Int, j::Int) = has_edge(net, j, i) ? 1.0 : 0.0
net = network(4; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 1); add_edge!(net, 3, 4)
change_stat(Mutuals(), net, 4, 3)                                    # 1.0
change_stat(Mutuals(), net, 1, 2)                                    # 1.0 — edge present, same value
change_stat_check(Mutuals(), net; n_tests=8, rng=Xoshiro(1))         # true
```
"""
change_stat

# =============================================================================
# Term Development Infrastructure
# =============================================================================

"""
    AbstractUserTerm <: AbstractERGMTerm

Base type for user-defined ERGM terms. A subtype is a full ERGM.jl term —
it can go into an `ERGMFormula`, be fitted with `fit_ergm`/`ergm` and
simulated from — once it has `name`, `compute` and `change_stat` methods
(the add-direction convention) and declares the traits that differ from the
defaults (`ERGM.required_vertex_attributes`, `ERGM.requires_directed`,
`ERGM.is_dyad_dependent`, `Networks.supports_missing`, …). The fallback
`ERGM.is_dyad_dependent(::AbstractUserTerm)` is `true`, the conservative
answer; ERGM's own `NodalTerm`/`DyadicTerm` supertypes imply `false`, so
subtype those only for covariate-only terms.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
import ERGM: name, compute, change_stat
struct EdgeCount <: AbstractUserTerm end
name(::EdgeCount) = "edgecount"
compute(::EdgeCount, net) = Float64(ne(net))
change_stat(::EdgeCount, net, i::Int, j::Int) = 1.0
ERGM.is_dyad_dependent(::EdgeCount) = false          # covariate-only: exact MPLE
EdgeCount() isa AbstractERGMTerm                     # true — a first-class ERGM term
net = network(6; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
validate_term(EdgeCount(), net; verbose=false, rng=Xoshiro(1))   # true
ERGMModel(ERGMFormula([EdgeCount()]), net) isa ERGMModel         # true
```
"""
abstract type AbstractUserTerm <: AbstractERGMTerm end

"""
    @ergm_term name body

Macro for defining custom ERGM terms. After evaluating `body`, it verifies
that the named type exists and that `compute`, `change_stat`, and `name`
methods are defined for it, warning about anything missing. (Numeric
consistency is checked separately with [`validate_term`](@ref) /
[`change_stat_check`](@ref), which need a term *instance* and a network.)

The generics must be imported by name before the macro is used (`import
ERGM: name, compute, change_stat`); with a bare `using`, the body defines
local functions of those names and the macro reports them as missing.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
import ERGM: name, compute, change_stat

@ergm_term ScaledEdges begin
    struct ScaledEdges <: AbstractUserTerm
        scale::Float64
    end

    name(t::ScaledEdges) = "scalededges.\$(t.scale)"

    compute(t::ScaledEdges, net) = t.scale * ne(net)

    # Add-direction: one more edge adds `scale`, whatever the dyad's state
    change_stat(t::ScaledEdges, net, i::Int, j::Int) = t.scale
end

net = network(6; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
compute(ScaledEdges(2.0), net)                                          # 4.0
change_stat(ScaledEdges(2.0), net, 1, 2)                                # 2.0
validate_term(ScaledEdges(2.0), net; verbose=false, rng=Xoshiro(1))     # true
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

# The most-documented mistake: `using ERGM` (no `import ERGM: …`) followed by
# `compute(::MyTerm, net) = …` silently defines a LOCAL `compute` that the
# shared generic never dispatches to. From inside the harness the symptom is
# indistinguishable from "no method at all", so every "no method" report
# carries the import line the author is missing.
const _IMPORT_HINT = "if you defined one after a bare `using`, it is a local " *
                     "function ERGM.jl never calls: add " *
                     "`import ERGM: name, compute, change_stat` before the " *
                     "definitions"
_no_method_hint(f::Symbol, ::Type{T}) where T =
    " — no $f method for $(nameof(T)) reaches the shared generic; $_IMPORT_HINT"

# Interface checks run by @ergm_term after the term is defined
function _check_term_definition(::Type{T}) where T
    if !(T <: AbstractERGMTerm)
        @warn "$(nameof(T)) does not subtype AbstractERGMTerm/AbstractUserTerm; " *
              "it will not work with ERGM.jl's TermSet"
        return nothing
    end
    _only_fallback(compute, Tuple{T, Any}) &&
        @warn "$(nameof(T)) has no compute(term, net) method" *
              _no_method_hint(:compute, T)
    _only_fallback(change_stat, Tuple{T, Any, Int, Int}) &&
        @warn "$(nameof(T)) has no change_stat(term, net, i, j) method" *
              _no_method_hint(:change_stat, T)
    _only_fallback(name, Tuple{T}) &&
        @warn "$(nameof(T)) has no name(term) method; ERGM will use the " *
              "default \"$(string(nameof(T)))\"-style name" *
              _no_method_hint(:name, T)
    return nothing
end
_check_term_definition(@nospecialize(x)) =
    @warn "@ergm_term body did not define a type (got $(typeof(x)))"

# =============================================================================
# Validation and Testing
# =============================================================================
#
# Randomness contract (root CLAUDE.md, panel 2026-09 item 6): every harness
# entry point takes `rng::AbstractRNG` (default `Random.default_rng()`) and
# draws nothing from anywhere else — no bare `rand`, no rng-less `randperm`,
# no `Random.seed!` (the test suite greps this file for them). Two calls with
# equal rngs visit the same dyads and return the same verdicts, at any thread
# count.
#
# Replayability: each entry point snapshots the rng's state on entry
# (`_HarnessContext`) and every failure message ends with a literal
# `rng=Xoshiro(0x…, 0x…, 0x…, 0x…)` that replays THAT call exactly — the
# harness draws straight from `rng`, so passing the snapshot back reproduces
# the same dyad sequence. (`TaskLocalRNG`, the default, copies to a `Xoshiro`;
# a `MersenneTwister` has no short literal, so its failures carry no hint —
# the caller who passed it already holds the seed.)

"""
    _HarnessContext(entry::Symbol, rng)

What a failure message needs to make the run reproducible: the public entry
point whose replay reproduces it and a literal of the rng state on entry (or
`nothing` when the rng type has no short literal). Threaded through the
private checks as an argument — never a global.
"""
struct _HarnessContext
    entry::Symbol
    rng_literal::Union{Nothing, String}
end
_HarnessContext(entry::Symbol, rng::AbstractRNG) =
    _HarnessContext(entry, _rng_literal(rng))

_hex(x::UInt64) = "0x" * string(x; base=16, pad=16)
_rng_literal(rng::Xoshiro) =
    "Xoshiro(" * join((_hex(UInt64(s)) for s in (rng.s0, rng.s1, rng.s2, rng.s3)), ", ") * ")"
_rng_literal(rng::TaskLocalRNG) = _rng_literal(copy(rng))
_rng_literal(::AbstractRNG) = nothing

# The suffix appended to every failure `@warn`; empty when no literal exists
function _replay_hint(ctx::_HarnessContext)
    ctx.rng_literal === nothing && return ""
    return "; reproduce with rng=$(ctx.rng_literal) passed to $(ctx.entry)"
end

# `showerror` text of a caught exception (an `ArgumentError`'s message, not
# the `ArgumentError("…")` constructor form `string(e)` would print)
_error_text(e) = sprint(showerror, e)

# A network ERGM.jl cannot fit is refused UP FRONT by every entry point that
# takes one, with ERGM's own message — never validated dyad by dyad and then
# rejected at the ERGMModel step in words that blame the term. Three cases:
#
# - fewer than 2 vertices: there is no dyad to check, so there is nothing a
#   verdict could rest on (`benchmark_term` has always refused this);
# - two-mode: `ERGMModel` refuses one, and every harness draw/enumeration
#   ranges over ALL off-diagonal pairs — on a two-mode network the
#   within-mode pairs are structurally absent (`add_edge!` returns `false`),
#   so the brute-force reference would read 0 there and a CORRECT term would
#   be reported inconsistent;
# - a self-loop present: `ERGMModel` refuses it because the term statistics
#   would count the loop while the estimators and samplers never touch the
#   diagonal. The harness would happily validate `compute` counting the loop.
#
# The two-mode and self-loop refusals ARE ERGM's — its `public`
# `_refuse_two_mode` / `_refuse_self_loops`, the functions `ERGMModel` runs —
# prefixed with the entry point, so the wording cannot drift from ERGM's
# (it used to be mirrored here by hand).
function _refuse_unfittable(net::Network, entry::Symbol)
    n = Int(nv(net))
    n >= 2 || throw(ArgumentError(
        "$entry needs a network with at least 2 vertices to draw dyads from, " *
        "got $n"))
    try
        ERGM._refuse_two_mode(net)
        ERGM._refuse_self_loops(net)
    catch e
        e isa ArgumentError || rethrow()
        throw(ArgumentError(
            "$entry validates terms on networks ERGMModel accepts only, and " *
            "ERGMModel refuses this one: " * e.msg))
    end
    return net
end

# Every network-taking entry point needs at least one dyad to draw; a
# `n_tests` of 0 would mean "skip the only numeric check" and hand out a PASS
# verdict resting on nothing.
function _require_positive(entry::Symbol, kw::Symbol, value::Int)
    value >= 1 || throw(ArgumentError(
        "$entry: $kw must be at least 1, got $value (a run that checks no dyad " *
        "cannot return a verdict)"))
    return value
end

# A term that does not declare `supports_missing` is validated at the face
# value of every masked dyad: that is what its statistic computes, and
# ERGM.jl applies the caller's `missing=` policy at estimation time, not the
# harness. The missing-data contract asks that a face-value number never be
# handed out without saying so — this is the saying so (logged only when the
# entry point is verbose; the returned verdict carries no number).
function _masked_face_value_notice(term::AbstractERGMTerm, net::Network)
    supports_missing(term) && return nothing
    k = n_missing_dyads(net)
    k == 0 && return nothing
    @info "net has $k masked dyad$(k == 1 ? "" : "s"); the term declares " *
          "supports_missing = false, so they are validated at their face value " *
          "(ERGM.jl applies its missing= policy at estimation time)"
    return nothing
end

"""
    validate_term(term::AbstractERGMTerm, net::Network; verbose=true, traits=true,
                  n_tests=10, rng=Random.default_rng()) -> Bool

Validate that a term is correctly implemented.

Checks the **method interface**:
- `compute()` returns a Real
- `change_stat()` returns a Real on `n_tests` random dyads
- `name()` returns a non-empty string
- Change statistics are consistent with compute differences (add-direction
  convention, [`change_stat_check`](@ref) on `n_tests` random dyads)

and, unless `traits=false`, the term's **trait declarations** — see
[`validate_traits`](@ref) for exactly what is asserted about
`ERGM.required_vertex_attributes`, `ERGM.required_edge_attributes`,
`ERGM.requires_directed`/`requires_undirected`, `ERGM.is_dyad_dependent` and
`Networks.supports_missing`. The trait checks also assert that the term is
accepted by `ERGMModel` construction on `net`, so a term that passes
`validate_term` is a term ERGM.jl will fit.

`net` must carry every vertex/edge attribute the term declares (validate
against a network the term is meant for), with a value on every vertex —
ERGM.jl refuses partial attributes at model construction, and so does this —
and with values that *vary* across vertices (a constant declared attribute
cannot be shown to be read, and the run fails). A network ERGM.jl cannot fit
is refused up front with an `ArgumentError` rather than validated: fewer
than 2 vertices (no dyad to check), a two-mode network, or one containing a
self-loop. A verdict is never returned without evidence: `n_tests` must be
at least 1.

# Keywords
- `verbose=true`: log each check (`@info`) and each failure (`@warn`).
- `traits=true`: run [`validate_traits`](@ref) after the interface checks.
- `n_tests=10`: random dyads drawn for the `change_stat` checks (at least 1;
  capped at the number of dyads in `net`).
- `rng=Random.default_rng()`: the source of **every** random draw. Two calls
  with equal rngs visit the same dyads. Every failure message ends with
  `reproduce with rng=Xoshiro(0x…, 0x…, 0x…, 0x…)`, the rng state on entry —
  pass it back to replay the run.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(10; directed=true)
for (i, j) in [(1, 2), (2, 3), (3, 1), (1, 4), (4, 5)]
    add_edge!(net, i, j)
end
validate_term(ExampleTerm(), net; verbose=false, rng=Xoshiro(1))   # true
```
"""
function validate_term(term::AbstractERGMTerm, net::Network; verbose::Bool=true,
                       traits::Bool=true, n_tests::Int=10,
                       rng::AbstractRNG=Random.default_rng())
    _require_positive(:validate_term, :n_tests, n_tests)
    _refuse_unfittable(net, :validate_term)
    ctx = _HarnessContext(:validate_term, rng)
    hint = _replay_hint(ctx)
    valid = true
    T = typeof(term)
    # "No method" is what a bare-`using` definition looks like from here
    no_compute = _only_fallback(compute, Tuple{T, Any})
    no_change = _only_fallback(change_stat, Tuple{T, Any, Int, Int})
    no_name = _only_fallback(name, Tuple{T})

    # Check name
    try
        n = name(term)
        if isempty(n)
            verbose && @warn "Term name is empty"
            valid = false
        elseif no_name
            # ERGM's fallback name is legitimate, but an author who just
            # wrote a name method wants to know it is not the one in use
            verbose && @info "✓ name() returns: $n (ERGM's fallback" *
                             _no_method_hint(:name, T) * ")"
        else
            verbose && @info "✓ name() returns: $n"
        end
    catch e
        verbose && @warn "name() failed: $(_error_text(e))" *
                         (no_name ? _no_method_hint(:name, T) : "")
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
        verbose && @warn "compute() failed: $(_error_text(e))" *
                         (no_compute ? _no_method_hint(:compute, T) : "")
        valid = false
    end
    verbose && _masked_face_value_notice(term, net)

    # Validate change statistics (`_refuse_unfittable` guarantees a dyad to
    # draw, `_require_positive` that at least one is drawn)
    n_vertices = Int(nv(net))
    n_draws = min(n_tests, _n_dyads(net))
    verbose && @info "Testing change_stat() with $n_draws random dyads..."

    for _ in 1:n_draws
        i, j = _random_dyad(rng, n_vertices)

        try
            delta = change_stat(term, net, i, j)
            if !isa(delta, Real)
                verbose && @warn "change_stat() should return a Real at ($i,$j), " *
                                 "got $(typeof(delta))$hint"
                valid = false
                break
            end
        catch e
            verbose && @warn "change_stat($i, $j) failed: $(_error_text(e))" *
                             (no_change ? _no_method_hint(:change_stat, T) : "") *
                             hint
            valid = false
            break
        end
    end

    if valid
        verbose && @info "✓ change_stat() returns valid values"
    end

    # Consistency check
    if valid
        verbose && @info "Running consistency check..."
        consistent = change_stat_check(term, net; n_tests=n_draws, verbose=false,
                                       rng=rng)
        if !consistent
            verbose && @warn "change_stat() values inconsistent with compute() " *
                             "differences (run change_stat_check with verbose=true " *
                             "for the dyad)$hint"
            valid = false
        else
            verbose && @info "✓ change_stat() consistent with compute()"
        end
    end

    # Trait declarations (the public protocol in ERGM.jl `src/terms/traits.jl`)
    if valid && traits
        valid &= validate_traits(term, net; verbose=verbose, rng=rng)
    end

    return valid
end

# Number of off-diagonal dyads of `net` (ordered on a directed network)
function _n_dyads(net::Network)
    n = Int(nv(net))
    return is_directed(net) ? n * (n - 1) : n * (n - 1) ÷ 2
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

# Vertices of `net` with no value for vertex attribute `attr` — absent from
# the attribute Dict or stored as `missing`/`nothing` (an NA on the way in
# from R). ERGM.jl refuses such an attribute at model construction (statnet
# refuses NA), so a declared attribute must be complete.
function _vertices_without_value(net::Network, attr::Symbol)
    raw = get_vertex_attribute(net, attr)
    return [v for v in 1:Int(nv(net))
            if (val = get(raw, v, nothing); val === nothing || val === missing)]
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
    validate_traits(term::AbstractERGMTerm, net::Network; verbose=true,
                    rng=Random.default_rng()) -> Bool

Exercise a term's declarations under ERGM.jl's public term-trait protocol
(`ERGM.required_vertex_attributes`, `ERGM.required_edge_attributes`,
`ERGM.requires_directed`, `ERGM.requires_undirected`,
`ERGM.is_dyad_dependent`, `Networks.supports_missing`) against `net`.

Checked, in order:

1. **Attributes are the ones it reads.** Every attribute the term declares
   must exist on `net` — a declared *vertex* attribute with a value on
   **every vertex**, since ERGM.jl refuses a partial attribute at model
   construction the way statnet refuses NA — and must actually move the
   term's statistic when its values are perturbed (otherwise the declaration
   is a warning: an over-declaration needlessly narrows the networks the term
   can be fitted to). Conversely, a term whose statistic moves when an
   **undeclared vertex attribute** is perturbed fails: it reads an attribute
   ERGM.jl will not validate, so on a network lacking that attribute it would
   silently become an all-zero design column instead of an error. (An
   undeclared *edge* attribute is only warned about — terms like
   [`WeightedEdges`](@ref) read one but fall back to a default, so they
   genuinely do not require it.) An attribute that is **constant** on `net`
   cannot be perturbed, so whether the term reads it is undecidable there:
   the check is reported as not run (never as ✓), and for a *declared*
   attribute the run fails — validate on a network where the attribute
   varies.
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

Finally, `ERGMModel` construction with the term on `net` must succeed
(`ERGM._validate_formula` is what runs there; its message is reported
verbatim on failure).

# Keywords
- `verbose=true`: log each check (`@info`) and each failure (`@warn`).
- `rng=Random.default_rng()`: the source of every random draw (the source
  dyads of check 3 and the masked dyads of check 4). Every failure message
  ends with `reproduce with rng=Xoshiro(0x…, …)`, the rng state on entry.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(6; directed=true)
set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:6))
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
validate_traits(InteractionTerm(:a, :a), net; verbose=false, rng=Xoshiro(1))  # true
```
"""
function validate_traits(term::AbstractERGMTerm, net::Network; verbose::Bool=true,
                         rng::AbstractRNG=Random.default_rng())
    _refuse_unfittable(net, :validate_traits)
    ctx = _HarnessContext(:validate_traits, rng)
    hint = _replay_hint(ctx)
    valid = true
    n = Int(nv(net))

    # --- 1. Attribute declarations ------------------------------------------
    undecided = Symbol[]      # attributes constant on `net`: dependence untestable
    for (kind, declared, present) in
            ((:vertex, ERGM.required_vertex_attributes(term), list_vertex_attributes(net)),
             (:edge, ERGM.required_edge_attributes(term), list_edge_attributes(net)))
        for attr in declared
            if !(attr in present)
                verbose && @warn "term declares required $kind attribute :$attr, " *
                                 "which the validation network does not have; " *
                                 "ERGMModel construction would reject it"
                valid = false
            elseif kind === :vertex
                absent = _vertices_without_value(net, attr)
                if !isempty(absent)
                    verbose && @warn "term declares required vertex attribute :$attr, " *
                                     "which must have a value on every vertex, but " *
                                     "$(length(absent)) of $n vertices have none " *
                                     "(vertices $(join(absent, ", "))); ERGMModel " *
                                     "construction refuses partial attributes as " *
                                     "statnet refuses NA — set a value on every vertex"
                    valid = false
                end
            end
        end
        for attr in present
            reads = _reads_attribute(term, net, attr, kind)
            if reads === missing
                # Constant on `net`: perturbing it changes nothing, so
                # whether the term reads it cannot be decided here. Nothing
                # was tested, so nothing may be reported as ✓ — and for a
                # DECLARED attribute the network is not one "the term is
                # meant for" (the docstring's requirement), so the run fails.
                push!(undecided, attr)
                if attr in declared
                    verbose && @warn "term declares required $kind attribute :$attr, " *
                                     "which is constant on this network, so whether " *
                                     "the term reads it is undecidable; validate on " *
                                     "a network where it varies"
                    valid = false
                else
                    verbose && @warn "$kind attribute :$attr is constant on this " *
                                     "network, so whether the term reads it is " *
                                     "undecidable; validate on a network where it " *
                                     "varies"
                end
                continue
            end
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
    if valid && verbose && isempty(undecided)
        @info "✓ declared attributes match the ones the term reads"
    elseif valid && verbose
        @info "declared attributes match the ones the term reads, except that " *
              "$(join(":" .* string.(undecided), ", ")) could not be tested " *
              "(constant on this network)"
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
    if !ERGM.is_dyad_dependent(term)
        work = deepcopy(net)
        sources = unique(_random_dyad(rng, n) for _ in 1:4)
        toggled_any = false     # a ✓ is printed only if some dyad was toggled
        for (i, j) in sources
            baseline = change_stat(term, work, i, j)
            for (k, l) in _other_dyads(rng, work, i, j)
                toggled_any = true
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
                                     "toggled: it IS dyad-dependent$hint"
                    valid = false
                    break
                end
            end
            valid || break
        end
        if valid && verbose && toggled_any
            @info "✓ is_dyad_dependent = false holds under toggling of other dyads"
        elseif valid && verbose
            # a 2-vertex undirected network has no other dyad to toggle
            @info "is_dyad_dependent = false could not be tested: the network " *
                  "has no other dyad to toggle"
        end
    end

    # --- 4. Missing-data declaration ----------------------------------------
    if supports_missing(term)
        work = deepcopy(net)
        for _ in 1:5
            i, j = _random_dyad(rng, n)
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
                                 "unobserved dyads at face value$hint"
                valid = false
                break
            end
        end
        valid && verbose && @info "✓ supports_missing = true holds: masked dyads do " *
                                  "not enter the statistic"
    end

    # --- 5. ERGM model construction accepts the term -------------------------
    # `ERGMModel` runs ERGM._validate_formula on the declarations: attribute
    # existence AND completeness (NA-refusal), direction, covariate sizes.
    # (The network-only refusals — two-mode, self-loops — were raised up
    # front by `_refuse_unfittable`, so a rejection here is about the term's
    # declarations on this network; the wording still leaves both open.)
    try
        ERGMModel(ERGMFormula([Edges(), term]), net)
        verbose && @info "✓ accepted by ERGMModel construction"
    catch e
        verbose && @warn "ERGMModel construction rejected the term or the " *
                         "network: $(_error_text(e))"
        valid = false
    end

    return valid
end

# Every dyad of `net` other than `(i, j)` itself (and, on an undirected
# network, its mirror). Capped at `cap` dyads chosen with `rng` so the check
# stays cheap on large networks.
function _other_dyads(rng::AbstractRNG, net::Network, i::Int, j::Int; cap::Int=250)
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
    length(out) > cap && (out = out[randperm(rng, length(out))[1:cap]])
    return out
end

# Uniform draw from the off-diagonal dyads (i ≠ j), without discarding draws
function _random_dyad(rng::AbstractRNG, n::Int)
    i = rand(rng, 1:n)
    j = rand(rng, 1:(n - 1))
    j >= i && (j += 1)
    return i, j
end

# Edge attributes are treated as exogenous dyadic data by the harness:
# Networks.rem_edge! deletes them, so they are snapshotted before a toggle
# and restored afterwards. Without this, any term reading edge attributes
# would spuriously fail the brute-force comparison. Both go through the
# public attribute API (`require_edge=false` stores a value ahead of the
# edge, which is exactly what restoring on a removed dyad needs).
function _edge_attr_snapshot(net::Network, i::Int, j::Int)
    key = is_directed(net) ? (i, j) : minmax(i, j)
    return [attr => attrs[key] for attr in list_edge_attributes(net)
            for attrs in (get_edge_attribute(net, attr),) if haskey(attrs, key)]
end

function _edge_attr_restore!(net::Network, i::Int, j::Int, saved)
    for (attr, val) in saved
        set_edge_attribute!(net, attr, i, j, val; require_edge=false)
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
# dyad's current state. `ctx` supplies the replay hint for failure messages.
function _dyad_consistent(term::AbstractERGMTerm, test_net::Network,
                          i::Int, j::Int, tol::Float64, verbose::Bool,
                          ctx::_HarnessContext)
    predicted = change_stat(term, test_net, i, j)
    expected = _brute_change_stat(term, test_net, i, j)

    if abs(predicted - expected) > tol
        verbose && @warn "Inconsistency at dyad ($i,$j): " *
                         "change_stat=$predicted, brute-force=$expected" *
                         _replay_hint(ctx)
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
                         "of whether the edge exists" * _replay_hint(ctx)
        return false
    end

    return true
end

"""
    change_stat_check(term::AbstractERGMTerm, net::Network; n_tests=10, verbose=true,
                      tol=1e-10, rng=Random.default_rng()) -> Bool

Verify that `change_stat` returns the **add-direction** change statistic
required by ERGM.jl: `compute` with edge (i,j) present minus `compute` with
it absent, independent of the dyad's current state. `n_tests` random dyads
(at least 1; drawn from `rng`, never `i == j`) are toggled on a copy of
`net` and restored; each is checked twice — as-is, and with the dyad's state
flipped — so a toggle-direction term fails the second check. A network with
fewer than 2 vertices, a two-mode one, or one containing a self-loop is
refused with an `ArgumentError` (ERGM.jl cannot fit it; see
[`validate_term`](@ref)).

Every failure message (`verbose=true`) names the dyad and ends with
`reproduce with rng=Xoshiro(0x…, …)`, the rng state on entry; pass it back
to replay the same dyad sequence.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(8; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3); add_edge!(net, 3, 1)
change_stat_check(ExampleTerm(), net; n_tests=5, rng=Xoshiro(1))          # true
change_stat_check(ExampleTerm(), net; n_tests=5, rng=Xoshiro(1)) ==
    change_stat_check(ExampleTerm(), net; n_tests=5, rng=Xoshiro(1))      # true
```
"""
function change_stat_check(term::AbstractERGMTerm, net::Network;
                           n_tests::Int=10, verbose::Bool=true, tol::Float64=1e-10,
                           rng::AbstractRNG=Random.default_rng())
    _require_positive(:change_stat_check, :n_tests, n_tests)
    _refuse_unfittable(net, :change_stat_check)
    verbose && _masked_face_value_notice(term, net)
    ctx = _HarnessContext(:change_stat_check, rng)
    test_net = deepcopy(net)
    n = Int(nv(test_net))
    all_passed = true

    for _ in 1:n_tests
        i, j = _random_dyad(rng, n)
        all_passed &= _dyad_consistent(term, test_net, i, j, tol, verbose, ctx)
    end

    return all_passed
end

"""
    consistency_check(term::AbstractERGMTerm, net::Network; exhaustive=false,
                      tol=1e-10, verbose=false, rng=Random.default_rng()) -> Bool

Full consistency check between `compute` and the add-direction
`change_stat` (see [`change_stat_check`](@ref)); stops at the first
inconsistent dyad. With `exhaustive=true` every dyad of `net` is checked
once — each ordered pair on a directed network, each unordered pair
(`i < j`) on an undirected one (can be slow for large networks); otherwise
up to `min(100, number of dyads)` distinct random dyads drawn from `rng`
(never `i == j`; canonicalised to `(min, max)` on an undirected network).
`verbose=true` reports the failing dyad, with the `rng=Xoshiro(0x…, …)`
replay literal. A network with fewer than 2 vertices, a two-mode one, or
one containing a self-loop is refused with an `ArgumentError`.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(6; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
consistency_check(TemplateTerm(2.0), net; exhaustive=true)           # true
consistency_check(TemplateTerm(2.0), net; rng=Xoshiro(3))            # true
```
"""
function consistency_check(term::AbstractERGMTerm, net::Network;
                           exhaustive::Bool=false, tol::Float64=1e-10,
                           verbose::Bool=false,
                           rng::AbstractRNG=Random.default_rng())
    _refuse_unfittable(net, :consistency_check)
    verbose && _masked_face_value_notice(term, net)
    ctx = _HarnessContext(:consistency_check, rng)
    n = Int(nv(net))
    test_net = deepcopy(net)
    directed = is_directed(net)

    # On an undirected network (i, j) and (j, i) are the same dyad: enumerate
    # each unordered dyad once (the `j < i` skip `_term_fingerprint` uses)
    # and size the random budget by the dyads that exist (`_n_dyads`), so
    # the check costs one brute-force comparison per dyad, not two
    dyads_to_check = if exhaustive
        [(i, j) for i in 1:n for j in 1:n if i != j && (directed || i < j)]
    else
        unique(directed ? _random_dyad(rng, n) : minmax(_random_dyad(rng, n)...)
               for _ in 1:min(100, _n_dyads(net)))
    end

    for (i, j) in dyads_to_check
        _dyad_consistent(term, test_net, i, j, tol, verbose, ctx) || return false
    end

    return true
end

"""
    test_term(term::AbstractERGMTerm; n_vertices=20, density=0.1, n_tests=100,
              directed=true, vertex_attributes=Dict{Symbol,Any}(),
              rng=Random.default_rng()) -> Bool

Comprehensive test suite for an ERGM term: builds a random network of
`n_vertices` vertices at edge `density` (directed unless `directed=false`;
edges drawn from `rng`), gives it the vertex attributes in
`vertex_attributes`, runs [`validate_term`](@ref) on it with `n_tests`
random dyads and the same `rng`, then checks `compute` on an empty and on a
complete network of the same kind (same directedness and attributes). Prints
a report (including the `rng=Xoshiro(0x…, …)` literal that replays the whole
run) and returns whether everything passed.

The generated networks carry **no** vertex attributes unless you ask, and
are directed unless you ask: a term that declares
`ERGM.required_vertex_attributes` needs `vertex_attributes`, one that
declares `ERGM.requires_undirected` needs `directed=false` — otherwise
`validate_term` (correctly) reports the mismatch and the run FAILS.

# Keywords
- `n_vertices=20`, `density=0.1`: size (at least 2) and edge density (in
  `[0, 1]`) of the random network; anything else is an `ArgumentError`, not
  a report on a network with no dyad to check.
- `n_tests=100`: random dyads for the `change_stat` checks (at least 1;
  capped at the number of dyads).
- `directed=true`: build a directed (`true`) or undirected (`false`) network.
- `vertex_attributes=Dict{Symbol,Any}()`: `attr => spec` pairs set on every
  generated network; `spec` is either a function of the vertex index
  (`v -> value`) or a vector with one value per vertex (length ≥ `n_vertices`).
- `rng=Random.default_rng()`: the source of every random draw — the network
  and the dyads.

# Example
```julia
using ERGM, ERGMUserterms, Random
test_term(ExampleTerm(); n_vertices=8, density=0.2, n_tests=5, rng=Xoshiro(2))   # true
# A term declaring vertex attributes needs the generated network to carry them
test_term(InteractionTerm(:a, :b); n_vertices=8, n_tests=5, rng=Xoshiro(2),
          vertex_attributes=Dict(:a => (v -> Float64(v)), :b => (v -> Float64(9 - v))))   # true
# A term defined on undirected networks only (ERGM.requires_undirected) is
# tested on one
test_term(TemplateTerm(2.0); n_vertices=8, n_tests=5, directed=false, rng=Xoshiro(2))   # true
```
"""
function test_term(term::AbstractERGMTerm;
                   n_vertices::Int=20,
                   density::Float64=0.1,
                   n_tests::Int=100,
                   directed::Bool=true,
                   vertex_attributes::AbstractDict{Symbol}=Dict{Symbol, Any}(),
                   rng::AbstractRNG=Random.default_rng())
    # Refuse a run that could only report a verdict resting on nothing:
    # `validate_term` would throw on these anyway, but after the header
    n_vertices >= 2 || throw(ArgumentError(
        "test_term needs a network with at least 2 vertices to draw dyads " *
        "from, got n_vertices=$n_vertices"))
    _require_positive(:test_term, :n_tests, n_tests)
    ctx = _HarnessContext(:test_term, rng)

    # Create test network (`_random_network` refuses a density outside [0, 1])
    net = _random_network(rng, n_vertices, density; directed=directed)
    _apply_vertex_attributes!(net, vertex_attributes)

    println("Testing term: $(name(term))")
    ctx.rng_literal === nothing ||
        println("(replay with rng=$(ctx.rng_literal))")
    println("=" ^ 50)

    # Run validation
    valid = validate_term(term, net; verbose=true, n_tests=n_tests, rng=rng)

    # Additional tests
    if valid
        println("\nAdditional tests:")

        # Test on empty network
        empty_net = network(n_vertices; directed=directed)
        _apply_vertex_attributes!(empty_net, vertex_attributes)
        try
            stat_empty = compute(term, empty_net)
            println("✓ Works on empty network: $stat_empty")
        catch e
            println("✗ Failed on empty network: $(_error_text(e))")
            valid = false
        end

        # Test on complete network
        m = min(10, n_vertices)
        complete_net = network(m; directed=directed)
        _apply_vertex_attributes!(complete_net, vertex_attributes)
        for i in 1:m, j in 1:m
            i == j && continue
            (!directed && j < i) && continue
            add_edge!(complete_net, i, j)
        end
        try
            stat_complete = compute(term, complete_net)
            println("✓ Works on complete network: $stat_complete")
        catch e
            println("✗ Failed on complete network: $(_error_text(e))")
            valid = false
        end
    end

    println("=" ^ 50)
    println(valid ? "All tests PASSED" : "Some tests FAILED")

    return valid
end

# A random network on `n` vertices with about `density` of its dyads tied,
# every draw from `rng` (shared by test_term and profile_term). `density`
# is a proportion of the dyads: anything outside [0, 1] is refused rather
# than silently clamped by `round`/re-drawn dyads.
function _random_network(rng::AbstractRNG, n::Int, density::Float64;
                         directed::Bool=true)
    0.0 <= density <= 1.0 || throw(ArgumentError(
        "density must lie in [0, 1] (a proportion of the network's dyads), " *
        "got $density"))
    net = network(n; directed=directed)
    n < 2 && return net
    n_edges = round(Int, density * (directed ? n * (n - 1) : n * (n - 1) ÷ 2))
    for _ in 1:n_edges
        i, j = _random_dyad(rng, n)
        add_edge!(net, i, j)
    end
    return net
end

# Set the `attr => spec` vertex attributes of `test_term` on `net`: `spec` is
# a function of the vertex index or a vector indexed by it.
function _apply_vertex_attributes!(net::Network, specs::AbstractDict{Symbol})
    n = Int(nv(net))
    for (attr, spec) in specs
        if spec isa AbstractVector && length(spec) < n
            throw(ArgumentError("vertex_attributes[:$attr] has $(length(spec)) " *
                                "values for a network of $n vertices"))
        end
        set_vertex_attribute!(net, attr,
                              Dict(v => _attribute_value(spec, v) for v in 1:n))
    end
    return net
end
_attribute_value(spec::AbstractVector, v::Int) = spec[v]
_attribute_value(spec, v::Int) = spec(v)

"""
    benchmark_term(term::AbstractERGMTerm, net::Network; n_iter=1000,
                   rng=Random.default_rng()) -> NamedTuple

Benchmark `compute()` and `change_stat()` performance: `n_iter` timed calls
of each (at least 1), the `change_stat` calls on random dyads drawn from
`rng`. Returns `(compute_mean, compute_std, change_stat_mean,
change_stat_std, speedup)` in seconds; `speedup = compute_mean /
change_stat_mean`. A network with fewer than 2 vertices, a two-mode one, or
one containing a self-loop is refused with an `ArgumentError`, as by every
validator.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(10; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
b = benchmark_term(ExampleTerm(), net; n_iter=100, rng=Xoshiro(1))
b.speedup > 0     # true
```
"""
function benchmark_term(term::AbstractERGMTerm, net::Network; n_iter::Int=1000,
                        rng::AbstractRNG=Random.default_rng())
    _require_positive(:benchmark_term, :n_iter, n_iter)
    _refuse_unfittable(net, :benchmark_term)     # < 2 vertices, two-mode, loops
    n = Int(nv(net))

    # Benchmark compute
    compute_times = Float64[]
    for _ in 1:n_iter
        t = @elapsed compute(term, net)
        push!(compute_times, t)
    end

    # Benchmark change_stat
    change_times = Float64[]
    for _ in 1:n_iter
        i, j = _random_dyad(rng, n)
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
                 n_iter=200, directed=true, vertex_attributes=Dict{Symbol,Any}(),
                 rng=Random.default_rng()) -> Vector{NamedTuple}

Profile a term's `compute`/`change_stat` cost across network sizes.
Returns one [`benchmark_term`](@ref) result per size (with `n_vertices` and
the realized edge count `ne` added), useful for spotting change statistics
that do not scale (e.g. O(edges) instead of O(degree)). The random networks
and the benchmarked dyads are all drawn from `rng`.

The generated networks are built exactly as [`test_term`](@ref) builds
them: directed unless `directed=false`, and carrying the vertex attributes
in `vertex_attributes` (`attr => v -> value` or `attr => vector`). Without
them a term declaring `ERGM.requires_undirected` cannot be profiled at all,
and an attribute-declaring term is timed on its attribute-absent fallback
branch rather than its real cost.

# Keywords
- `sizes=[10, 20, 40]`: vertex counts (each at least 2) of the random networks.
- `density=0.1`: edge density of each random network, in `[0, 1]`.
- `n_iter=200`: timed calls per size (at least 1), passed to `benchmark_term`.
- `directed=true`: build directed (`true`) or undirected (`false`) networks.
- `vertex_attributes=Dict{Symbol,Any}()`: `attr => spec` pairs set on every
  generated network; `spec` is a function of the vertex index or a vector
  with one value per vertex (length ≥ the largest size).
- `rng=Random.default_rng()`: the source of every random draw.

# Example
```julia
using ERGM, ERGMUserterms, Random
prof = profile_term(ExampleTerm(); sizes=[5, 10], n_iter=5, rng=Xoshiro(1))
[r.n_vertices for r in prof]      # [5, 10]
# An attribute-declaring term is timed on networks that carry its attributes
prof = profile_term(InteractionTerm(:a, :b); sizes=[5, 10], n_iter=5, rng=Xoshiro(1),
                    vertex_attributes=Dict(:a => (v -> Float64(v)), :b => (v -> 1.0)))
length(prof)                      # 2
# ... and an undirected-only term on undirected networks
profile_term(TemplateTerm(1.0); sizes=[6], n_iter=5, directed=false, rng=Xoshiro(1))[1].ne <= 15   # true
```
"""
function profile_term(term::AbstractERGMTerm;
                      sizes::Vector{Int}=[10, 20, 40],
                      density::Float64=0.1,
                      n_iter::Int=200,
                      directed::Bool=true,
                      vertex_attributes::AbstractDict{Symbol}=Dict{Symbol, Any}(),
                      rng::AbstractRNG=Random.default_rng())
    results = NamedTuple[]
    for n in sizes
        net = _random_network(rng, n, density; directed=directed)
        _apply_vertex_attributes!(net, vertex_attributes)
        b = benchmark_term(term, net; n_iter=n_iter, rng=rng)
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

Identical to statnet's `nodecov("id")` with `id` = vertex index; the golden
fixture `test/fixtures/userterms_examples.toml` pins it against ergm on
every dyad. Covariate-only, so it declares `ERGM.is_dyad_dependent = false`.

# Example
```julia
using ERGM, ERGMUserterms, Networks
net = network(5; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 3, 4)
compute(ExampleTerm(), net)              # 10.0 — (1 + 2) + (3 + 4)
change_stat(ExampleTerm(), net, 2, 3)    # 5.0
change_stat(ExampleTerm(), net, 1, 2)    # 3.0 — edge present: same add-direction value
name(ExampleTerm())                      # "example"
```
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

A parameterized template term — `param` times the number of edges. Copy
and modify this for your own terms. `TemplateTerm(1.0)` is statnet's
`edges` (pinned by the golden fixture).

# Fields
- `param::T`: A parameter value
- `attr::Symbol`: An attribute name (optional; unused by the template's
  statistic — a placeholder for the attribute your term reads)

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(5; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
t = TemplateTerm(2.5)
name(t)                                  # "template.2.5"
compute(t, net)                          # 5.0 — 2 edges × 2.5
change_stat(t, net, 3, 4)                # 2.5
TemplateTerm(3; attr=:x).attr            # :x
validate_term(t, net; verbose=false, rng=Xoshiro(1))   # true
```
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

`compute` sums over the network's edges, reading each edge's weight with the
per-edge getter `get_edge_attribute(net, attr, i, j)` (which returns
`nothing` for an edge — or a network — without the attribute) and falling
back to `default`; that keeps it consistent with `change_stat` when edges
are added without weights (e.g. by the MH sampler). For the same reason the
term declares **no** required edge attribute: it is well defined on a
network that has never heard of `attr`. On an undirected network the getter
reads the canonical `(min, max)` key itself. Identical to statnet's
`edgecov(W)` with `W` the stored weights and `default` elsewhere.

The read is asserted `::Float64`, so `compute`/`change_stat` infer `Float64`
and allocate nothing per call (the attribute store is untyped, `Dict{…,Any}`;
without the assertion the term would infer `Any` and box on every MPLE row
and MH step of any model containing it).

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(4; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 3)
set_edge_attribute!(net, :weight, 1, 2, 3.0)
t = WeightedEdges()                      # attr = :weight, default = 1.0
compute(t, net)                          # 4.0 — 3.0 + the default 1.0
change_stat(t, net, 1, 2)                # 3.0
change_stat(t, net, 3, 4)                # 1.0 — a fresh edge gets the default
ERGM.required_edge_attributes(t)         # () — deliberately undeclared
validate_term(t, net; verbose=false, rng=Xoshiro(1))   # true
```
"""
struct WeightedEdges <: AbstractUserTerm
    attr::Symbol
    default::Float64
    WeightedEdges(attr::Symbol=:weight; default::Float64=1.0) = new(attr, default)
end

name(t::WeightedEdges) = "weightedges.$(t.attr)"

# The weight edge (i,j) carries: its stored attribute if present, otherwise
# the default a fresh edge would get. The per-edge getter canonicalises the
# key on undirected networks, returns `nothing` when the edge (or the whole
# attribute) has no value, and allocates nothing — unlike
# `get_edge_attribute(net, attr)`, whose `Dict{Tuple{Int,Int},Any}` is
# untyped (a bare `get` from it infers `Any`, which would poison ERGM's
# change-statistic tuple: `Tuple{Float64, Any}`, boxed on every MPLE row and
# MH step). The `::Float64` assertion closes inference.
function _edge_weight(t::WeightedEdges, net, i::Int, j::Int)
    w = get_edge_attribute(net, t.attr, i, j)
    return w === nothing ? t.default : Float64(w)::Float64
end

function compute(t::WeightedEdges, net)
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

"""
    DyadCovTerm <: AbstractUserTerm

Dyadic covariate term: sum over edges of `covariate[i, j]`. Demonstrates
using a matrix covariate. Identical to statnet's `edgecov(covariate)`; on an
undirected network the entry is read at the canonical `(min, max)` position,
so pass a symmetric matrix (or accept that only its upper triangle counts).
Entries outside the matrix contribute 0.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
distance = [0.0 1.0 2.0 3.0
            1.0 0.0 1.5 2.5
            2.0 1.5 0.0 1.0
            3.0 2.5 1.0 0.0]
net = network(4; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 3, 4)
t = DyadCovTerm(distance)
compute(t, net)                          # 2.0 — distance[1,2] + distance[3,4]
change_stat(t, net, 1, 3)                # 2.0
validate_term(t, net; verbose=false, rng=Xoshiro(1))   # true
```
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

Interaction between two vertex attributes: the sum over edges (i, j) of
`a_i b_j + a_j b_i`. Identical to statnet's `edgecov(M)` with
`M = a bᵀ + b aᵀ`. It **declares** both attributes
(`ERGM.required_vertex_attributes`), so `ERGMModel` construction refuses a
network on which either is absent or set on only some vertices — the
statnet NA rule — instead of fitting a silent all-zero column. Values are
read per vertex with `get_vertex_attribute(net, attr, v)` and asserted
`::Float64`, so `compute`/`change_stat` infer `Float64` and allocate nothing.

# Example
```julia
using ERGM, ERGMUserterms, Networks, Random
net = network(4; directed=true)
add_edge!(net, 1, 2)
set_vertex_attribute!(net, :age, Dict(1 => 25.0, 2 => 30.0, 3 => 40.0, 4 => 50.0))
set_vertex_attribute!(net, :income, Dict(1 => 5.0, 2 => 6.0, 3 => 7.0, 4 => 8.0))
t = InteractionTerm(:age, :income)
compute(t, net)                          # 300.0 — 25·6 + 30·5
ERGM.required_vertex_attributes(t)       # (:age, :income)
validate_term(t, net; verbose=false, rng=Xoshiro(1))   # true
try; ERGMModel(ERGMFormula([InteractionTerm(:age, :height)]), net); catch e; e isa ArgumentError; end   # true
```
"""
struct InteractionTerm <: AbstractUserTerm
    attr1::Symbol
    attr2::Symbol
end

name(t::InteractionTerm) = "interact.$(t.attr1).$(t.attr2)"

function compute(t::InteractionTerm, net)
    # Note: `get_vertex_attribute(net, attr, v)` returns `nothing` when the
    # attribute is absent or the vertex has no value, and `_vertex_value`
    # zero-fills. That zero-fill is unreachable through `ERGMModel`: the term
    # DECLARES both attributes (`ERGM.required_vertex_attributes` below), so
    # `ERGM._validate_formula` refuses a network on which either is absent or
    # set on only some vertices (statnet refuses NA) before any statistic is
    # computed. It exists only for raw `compute`/`change_stat` calls.
    total = 0.0
    for e in edges(net)
        total += _interaction(t, net, Int(src(e)), Int(dst(e)))
    end
    return total
end

# One vertex's value of `attr`, 0.0 when it has none. The per-vertex getter
# allocates nothing (unlike `get_vertex_attribute(net, attr)`, whose
# `Dict{Int,Any}` is untyped and would make the arithmetic below `Any * Any`
# — boxed on every call, and poisoning ERGM's change-statistic tuple); the
# `::Float64` assertion closes inference.
function _vertex_value(net, attr::Symbol, v::Int)
    x = get_vertex_attribute(net, attr, v)
    return x === nothing ? 0.0 : Float64(x)::Float64
end

# The per-edge contribution of (i, j): a_i b_j + a_j b_i
_interaction(t::InteractionTerm, net, i::Int, j::Int) =
    _vertex_value(net, t.attr1, i) * _vertex_value(net, t.attr2, j) +
    _vertex_value(net, t.attr1, j) * _vertex_value(net, t.attr2, i)

function change_stat(t::InteractionTerm, net, i::Int, j::Int)
    # Add-direction change: the per-edge contribution of (i,j)
    return _interaction(t, net, i, j)
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

# Example
```julia
using ERGM, ERGMUserterms
term_signature(TemplateTerm(2.5))        # "TemplateTerm(param::Float64, attr::Symbol)"
term_signature(ExampleTerm())            # "ExampleTerm"
```
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

Generate a Markdown documentation stub for a term: heading, signature
([`term_signature`](@ref)), `name`, one line per field with its current
value, and the three interface methods.

# Example
```julia
using ERGM, ERGMUserterms
doc = term_documentation(TemplateTerm(2.5))
occursin("Name: template.2.5", doc)                     # true
occursin("`param::Float64` = 2.5", doc)                 # true
```
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
