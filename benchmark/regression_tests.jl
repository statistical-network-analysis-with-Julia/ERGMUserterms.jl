#!/usr/bin/env julia
# benchmark/regression_tests.jl — allocation-regression assertions for the
# bundled ERGMUserterms.jl terms. Standalone; run with
#     julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
#     julia --project=benchmark benchmark/regression_tests.jl
#
# The same pins as test/runtests.jl's "Allocation regressions" testset, kept
# here so the site's tools/run_benchmarks.jl reports a row for this package
# and CI runs the benchmark environment itself (panel 2026-09, item 7).
#
# Every bundled term and the template's `ReciprocatedHomophily` is
# allocation-free per dyad and per full statistic, and `change_stat`/`compute`
# infer `Float64`. The attribute-reading terms get there by reading per
# element (`get_vertex_attribute(net, attr, v)` / `get_edge_attribute(net,
# attr, i, j)`, 0 B) and asserting `::Float64`; the whole-Dict getters
# return untyped `Dict{…,Any}`s (and allocate an empty default on every
# call), and a term reading them with a bare `get` infers `Any` — which is
# what round 1 shipped for `WeightedEdges`/`InteractionTerm` (96 B / 208 B
# per call, and ERGM's change-statistic tuple inferred `Tuple{Float64,
# Any}`). A typed snapshot at construction (`vertex_attribute_vector(net,
# attr, Float64)`) remains the faster pattern for a hot path. Any allocation
# or an `Any` return type is a regression.

using ERGM
using ERGMUserterms
using Networks
using Random
using Test
using Graphs: src, dst, edges

include(joinpath(@__DIR__, "..", "examples", "MyTermPackage", "src", "MyTermPackage.jl"))
using .MyTermPackage: ReciprocatedHomophily

function er_network(rng::AbstractRNG, n::Int, m::Int)
    net = network(n; directed=true)
    while ne(net) < m
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i == j && continue
        add_edge!(net, i, j)
    end
    set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:n))
    set_vertex_attribute!(net, :b, Dict(v => Float64(n + 1 - v) for v in 1:n))
    set_vertex_attribute!(net, :group, Dict(v => isodd(v) ? "a" : "b" for v in 1:n))
    for e in first(collect(edges(net)), 50)
        set_edge_attribute!(net, :weight, src(e), dst(e), 2.5)
    end
    return net
end

"Bytes allocated by `change_stat` on a pre-warmed call, worst over `dyads`."
function worst_change_alloc(term, net, dyads)
    worst = 0
    for (i, j) in dyads
        change_stat(term, net, i, j)                     # warm up / compile
        worst = max(worst, @allocated change_stat(term, net, i, j))
    end
    return worst
end

"Bytes allocated by `compute` on a pre-warmed call."
function compute_alloc(term, net)
    compute(term, net)
    return @allocated compute(term, net)
end

@testset "ERGMUserterms allocation regressions" begin
    n = 500
    net = er_network(Random.Xoshiro(20260912), n, 10 * n)
    rng = Random.Xoshiro(1)
    dyads = Tuple{Int, Int}[]
    while length(dyads) < 40
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i == j || push!(dyads, (i, j))
    end
    cov = [Float64(2i + j) for i in 1:n, j in 1:n]

    @testset "every bundled term is allocation-free and infers Float64" begin
        for term in (ExampleTerm(), TemplateTerm(1.5), WeightedEdges(), DyadCovTerm(cov),
                     InteractionTerm(:a, :b), ReciprocatedHomophily(:group))
            @test worst_change_alloc(term, net, dyads) == 0
            @test compute_alloc(term, net) == 0
            @test Base.return_types(change_stat, (typeof(term), typeof(net), Int, Int)) == [Float64]
            @test Base.return_types(compute, (typeof(term), typeof(net))) == [Float64]
        end
    end

    @testset "harness building blocks are capped" begin
        big = network(30; directed=true)
        for k in 1:29
            add_edge!(big, k, k + 1)
        end
        @test length(ERGMUserterms._other_dyads(Random.Xoshiro(1), big, 1, 2)) == 250
    end
end
