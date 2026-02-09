"""
    ERGMUserterms.jl - Custom ERGM Term Development

Provides templates, utilities, and validation tools for developing custom ERGM terms.
Includes example terms and a comprehensive testing framework.

Port of the R ergm.userterms package from the StatNet collection.
"""
module ERGMUserterms

using ERGM
using Graphs
using Network
using Random
using Statistics
using Test

# Term development
export @ergm_term, validate_term, test_term
export AbstractUserTerm

# Templates and examples
export ExampleTerm, TemplateTerm
export WeightedEdges, DyadCovTerm, InteractionTerm

# Testing utilities
export change_stat_check, consistency_check
export benchmark_term, profile_term

# Documentation helpers
export term_signature, term_documentation

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

Macro for defining custom ERGM terms with automatic validation.

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

        # Validate the term was defined correctly
        if isdefined(@__MODULE__, $(QuoteNode(termname)))
            @info "ERGM term $($(QuoteNode(termname))) defined"
        end
    end
end

# =============================================================================
# Validation and Testing
# =============================================================================

"""
    validate_term(term::AbstractERGMTerm, net::Network; verbose=true) -> Bool

Validate that a term is correctly implemented.

Checks:
- `compute()` returns a Real
- `change_stat()` returns correct values
- `name()` returns a non-empty string
- Change statistics are consistent with compute differences
"""
function validate_term(term::AbstractERGMTerm, net::Network; verbose::Bool=true)
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
            j = rand(1:n_vertices)
            i == j && continue

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

    return valid
end

"""
    change_stat_check(term::AbstractERGMTerm, net::Network; n_tests=10, verbose=true) -> Bool

Verify that change_stat() correctly predicts differences in compute().
Toggles edges and checks that change_stat matches compute(after) - compute(before).
"""
function change_stat_check(term::AbstractERGMTerm, net::Network;
                           n_tests::Int=10, verbose::Bool=true, tol::Float64=1e-10)
    test_net = deepcopy(net)
    n = nv(test_net)
    all_passed = true

    for _ in 1:n_tests
        i = rand(1:n)
        j = rand(1:n)
        i == j && continue

        # Compute statistic before
        stat_before = compute(term, test_net)

        # Get predicted change
        predicted_delta = change_stat(term, test_net, i, j)

        # Toggle edge
        if has_edge(test_net, i, j)
            rem_edge!(test_net, i, j)
        else
            add_edge!(test_net, i, j)
        end

        # Compute statistic after
        stat_after = compute(term, test_net)

        # Check consistency
        actual_delta = stat_after - stat_before
        if abs(predicted_delta - actual_delta) > tol
            if verbose
                @warn "Inconsistency at edge ($i,$j): " *
                      "predicted=$predicted_delta, actual=$actual_delta"
            end
            all_passed = false
        end

        # Toggle back
        if has_edge(test_net, i, j)
            rem_edge!(test_net, i, j)
        else
            add_edge!(test_net, i, j)
        end
    end

    return all_passed
end

"""
    consistency_check(term::AbstractERGMTerm, net::Network; exhaustive=false) -> Bool

Full consistency check between compute() and change_stat().
If exhaustive=true, checks all possible edges (can be slow for large networks).
"""
function consistency_check(term::AbstractERGMTerm, net::Network;
                           exhaustive::Bool=false, tol::Float64=1e-10)
    n = nv(net)
    test_net = deepcopy(net)

    edges_to_check = if exhaustive
        [(i, j) for i in 1:n for j in 1:n if i != j]
    else
        [(rand(1:n), rand(1:n)) for _ in 1:min(100, n*n)]
    end

    for (i, j) in edges_to_check
        i == j && continue

        stat_before = compute(term, test_net)
        delta = change_stat(term, test_net, i, j)

        if has_edge(test_net, i, j)
            rem_edge!(test_net, i, j)
        else
            add_edge!(test_net, i, j)
        end

        stat_after = compute(term, test_net)

        if abs(delta - (stat_after - stat_before)) > tol
            return false
        end

        # Restore
        if has_edge(test_net, i, j)
            rem_edge!(test_net, i, j)
        else
            add_edge!(test_net, i, j)
        end
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
        i, j = rand(1:n), rand(1:n)
        i == j && continue
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

# =============================================================================
# Example Terms and Templates
# =============================================================================

"""
    ExampleTerm <: AbstractUserTerm

Example custom term that counts edges weighted by vertex ID sum.
Use this as a template for creating your own terms.
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
    has_edge(net, i, j) ? -(i + j) : Float64(i + j)
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
    # Template: change is just the parameter (positive if adding, negative if removing)
    return has_edge(net, i, j) ? -t.param : t.param
end

"""
    WeightedEdges <: AbstractUserTerm

Sum of edge weights. Demonstrates accessing edge attributes.
"""
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
    i <= size(t.covariate, 1) && j <= size(t.covariate, 2) || return 0.0
    val = t.covariate[i, j]
    return has_edge(net, i, j) ? -val : val
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
