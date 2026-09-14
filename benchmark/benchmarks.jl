#!/usr/bin/env julia
# benchmark/benchmarks.jl — BenchmarkTools suite for ERGMUserterms.jl's
# bundled terms (the reference implementations third-party authors copy).
#
# One `@benchmarkable` per bundled term for `compute` (the full statistic,
# O(edges)) and for `change_stat` (the per-dyad add-direction change, which
# must be O(1) or O(degree) — never O(edges)), on sparse Erdős–Rényi networks
# with the SAME expected mean degree at n = 500 and n = 2000. The n-scaling of
# the measured per-dyad cost is asserted: an O(n)/O(edges) regression shows up
# as a ≈4× ratio, O(1)/O(degree) stays ≈1× (limit 3×).
#
# Defines the standard `SUITE::BenchmarkGroup`. Run standalone with
#     julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
#     julia --project=benchmark benchmark/benchmarks.jl
# which tunes + runs the suite, prints one tab-separated `BENCHJL` line per
# benchmark (consumed by the site repo's tools/run_benchmarks.jl), and exits
# non-zero if the scaling assertion fails. The allocation pins live in
# regression_tests.jl (and in test/runtests.jl, testset "Allocation
# regressions").

using BenchmarkTools
using ERGM
using ERGMUserterms
using Networks
using Random
using Graphs: src, dst, edges

# The package template's term (examples/MyTermPackage) is benchmarked with
# the bundled ones: it is the copyable skeleton
include(joinpath(@__DIR__, "..", "examples", "MyTermPackage", "src", "MyTermPackage.jl"))
using .MyTermPackage: ReciprocatedHomophily

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

const MEAN_DEGREE = 10          # identical at both sizes: cost must not grow
const N_SMALL = 500
const N_LARGE = 2000
const N_DYADS = 200             # dyads swept per benchmark evaluation
const SCALING_LIMIT = 3.0       # tolerated t(n=2000)/t(n=500) ratio
const N_WEIGHTED = 100          # edges carrying a stored :weight

"Sparse directed Erdős–Rényi network with expected mean degree `MEAN_DEGREE`,
carrying the attributes the bundled terms read (:a, :b, :group, :weight)."
function er_network(rng::AbstractRNG, n::Int)
    net = network(n; directed=true)
    m = MEAN_DEGREE * n
    while ne(net) < m
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i == j && continue
        add_edge!(net, i, j)
    end
    set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:n))
    set_vertex_attribute!(net, :b, Dict(v => Float64(n + 1 - v) for v in 1:n))
    set_vertex_attribute!(net, :group, Dict(v => isodd(v) ? "a" : "b" for v in 1:n))
    for e in first(collect(edges(net)), N_WEIGHTED)
        set_edge_attribute!(net, :weight, src(e), dst(e), 2.5)
    end
    return net
end

"Fixed sample of `N_DYADS` random dyads (i ≠ j) to sweep per evaluation."
function sample_dyads(rng::AbstractRNG, n::Int)
    dyads = Tuple{Int, Int}[]
    while length(dyads) < N_DYADS
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i == j || push!(dyads, (i, j))
    end
    return dyads
end

"Sum of add-direction change statistics over a fixed dyad sample."
function sweep_change_stat(term, net, dyads)
    s = 0.0
    for (i, j) in dyads
        s += change_stat(term, net, i, j)
    end
    return s
end

const NETS = Dict(n => er_network(Random.Xoshiro(n), n) for n in (N_SMALL, N_LARGE))
const DYADS = Dict(n => sample_dyads(Random.Xoshiro(n + 1), n) for n in (N_SMALL, N_LARGE))
const COVS = Dict(n => [Float64(2i + j) for i in 1:n, j in 1:n] for n in (N_SMALL, N_LARGE))

# Terms per network size (DyadCovTerm's matrix must match n)
terms_for(n) = [("example", ExampleTerm()),
                ("template", TemplateTerm(1.5)),
                ("weightededges", WeightedEdges()),
                ("dyadcov", DyadCovTerm(COVS[n])),
                ("interaction", InteractionTerm(:a, :b)),
                ("recip_homophily", ReciprocatedHomophily(:group))]
const TERM_LABELS = first.(terms_for(N_SMALL))

# ---------------------------------------------------------------------------
# Suite
# ---------------------------------------------------------------------------

const SUITE = BenchmarkGroup()

let g = addgroup!(SUITE, "compute"), h = addgroup!(SUITE, "change_stat")
    for n in (N_SMALL, N_LARGE), (label, term) in terms_for(n)
        g["$(label)_n$(n)"] = @benchmarkable compute($term, $(NETS[n]))
        h["$(label)_n$(n)"] = @benchmarkable sweep_change_stat($term, $(NETS[n]), $(DYADS[n]))
    end
end

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

function print_benchjl(results::BenchmarkGroup)
    for (path, trial) in BenchmarkTools.leaves(results)
        est = median(trial)
        println("BENCHJL\t", join(path, "/"), "\t",
                BenchmarkTools.time(est), "\t",
                BenchmarkTools.allocs(est), "\t",
                BenchmarkTools.memory(est))
    end
end

"Assert that per-dyad cost did not grow with n (O(1)/O(degree), not O(edges))."
function assert_scaling(results::BenchmarkGroup)
    ok = true
    for label in TERM_LABELS
        t_small = BenchmarkTools.time(median(results["change_stat"]["$(label)_n$(N_SMALL)"]))
        t_large = BenchmarkTools.time(median(results["change_stat"]["$(label)_n$(N_LARGE)"]))
        ratio = t_large / t_small
        println("SCALING\t", label, "\tn", N_LARGE, "/n", N_SMALL, "\t",
                round(ratio, digits=2))
        if ratio > SCALING_LIMIT
            println(stderr, "SCALING FAILURE: $label change statistic is ",
                    round(ratio, digits=2), "x slower at n=$(N_LARGE) than at ",
                    "n=$(N_SMALL) (same mean degree; limit $(SCALING_LIMIT)x). ",
                    "The per-dyad cost is no longer O(1)/O(degree).")
            ok = false
        end
    end
    return ok
end

function main()
    tune!(SUITE)
    results = run(SUITE; verbose=false, seconds=1)
    print_benchjl(results)
    assert_scaling(results) || exit(1)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
