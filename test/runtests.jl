using ERGMUserterms
using ERGM
using Network
using Graphs
using Random
using Test

import ERGMUserterms: name, compute, change_stat

# A random directed network fixture
function random_net(n, n_edges; directed=true, seed=42)
    rng = Random.Xoshiro(seed)
    net = network(n; directed=directed)
    while ne(net) < n_edges
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i != j && add_edge!(net, i, j)
    end
    return net
end

# A deliberately broken term using the toggle-direction convention; the
# harness must reject it
struct ToggleConventionTerm <: AbstractUserTerm end
name(::ToggleConventionTerm) = "toggle_broken"
compute(::ToggleConventionTerm, net) = Float64(ne(net))
change_stat(::ToggleConventionTerm, net, i::Int, j::Int) =
    has_edge(net, i, j) ? -1.0 : 1.0

# The doc walkthrough term (docs/src/getting_started.md), kept in sync so
# the documentation example is machine-verified
struct SharedNeighborTerm <: AbstractUserTerm end
name(::SharedNeighborTerm) = "shared_neighbors"
function compute(::SharedNeighborTerm, net)
    total = 0.0
    for e in edges(net)
        i, j = src(e), dst(e)
        for k in outneighbors(net, i)
            k != j && has_edge(net, j, k) && (total += 1.0)
        end
    end
    return total
end
function change_stat(::SharedNeighborTerm, net, i::Int, j::Int)
    delta = 0.0
    for k in outneighbors(net, i)
        k == j && continue
        has_edge(net, j, k) && (delta += 1.0)
    end
    for b in outneighbors(net, i)
        b == j && continue
        has_edge(net, b, j) && (delta += 1.0)
    end
    for a in inneighbors(net, i)
        a == j && continue
        has_edge(net, a, j) && (delta += 1.0)
    end
    return delta
end

# Terms for the @ergm_term definition checks (structs must be defined at
# top level)
struct CompleteTerm <: AbstractUserTerm end
name(::CompleteTerm) = "complete"
compute(::CompleteTerm, net) = Float64(ne(net))
change_stat(::CompleteTerm, net, i::Int, j::Int) = 1.0

struct IncompleteTerm <: AbstractUserTerm end

@testset "ERGMUserterms.jl" begin
    @testset "AbstractUserTerm hierarchy" begin
        @test AbstractUserTerm <: ERGM.AbstractERGMTerm
    end

    @testset "Term construction" begin
        @test ExampleTerm() isa AbstractUserTerm
        @test TemplateTerm(2.5).param == 2.5
        @test TemplateTerm(3; attr=:x).attr == :x
        @test WeightedEdges().attr == :weight
        @test WeightedEdges(:strength; default=0.0).default == 0.0
        @test DyadCovTerm(zeros(3, 3)) isa AbstractUserTerm
        @test InteractionTerm(:a, :b) isa AbstractUserTerm
    end

    @testset "name() extends ERGM.name" begin
        # Names must reach ERGM.jl's machinery (TermSet stores names via
        # ERGM.name)
        @test ERGM.name(ExampleTerm()) == "example"
        @test ERGM.name(TemplateTerm(2.5)) == "template.2.5"

        ts = TermSet([Edges(), ExampleTerm()])
        @test ts.names == ["edges", "example"]
    end

    @testset "Example terms satisfy the add-direction invariant" begin
        Random.seed!(7)
        for directed in (true, false)
            net = random_net(12, 20; directed=directed, seed=directed ? 1 : 2)
            set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:12))
            set_vertex_attribute!(net, :b, Dict(v => Float64(13 - v) for v in 1:12))
            for e in collect(edges(net))[1:5]
                set_edge_attribute!(net, :weight, src(e), dst(e), 2.5)
            end

            cov = [Float64(i * 2 + j) for i in 1:12, j in 1:12]  # asymmetric

            for term in (ExampleTerm(), TemplateTerm(1.5), WeightedEdges(),
                         DyadCovTerm(cov), InteractionTerm(:a, :b))
                @test consistency_check(term, net; exhaustive=true)
            end
        end
    end

    @testset "Harness rejects toggle-direction terms" begin
        net = random_net(10, 15)
        @test !change_stat_check(ToggleConventionTerm(), net;
                                 n_tests=20, verbose=false)
        @test !consistency_check(ToggleConventionTerm(), net; exhaustive=true)
    end

    @testset "Doc walkthrough term is correct" begin
        net = random_net(10, 18; seed=3)
        @test consistency_check(SharedNeighborTerm(), net; exhaustive=true)
    end

    @testset "validate_term" begin
        net = random_net(10, 15)
        @test validate_term(ExampleTerm(), net; verbose=false)
        @test validate_term(WeightedEdges(), net; verbose=false)
        @test !validate_term(ToggleConventionTerm(), net; verbose=false)
    end

    @testset "WeightedEdges semantics" begin
        net = network(4)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        set_edge_attribute!(net, :weight, 1, 2, 3.0)

        term = WeightedEdges()
        # Stored weight + default for the unweighted edge
        @test compute(term, net) == 3.0 + 1.0
        @test change_stat(term, net, 1, 2) == 3.0
        @test change_stat(term, net, 3, 4) == 1.0   # fresh edge gets default

        # No weight attribute at all: every edge counts at the default
        bare = network(3)
        add_edge!(bare, 1, 2)
        @test compute(term, bare) == 1.0
        @test consistency_check(term, bare; exhaustive=true)
    end

    @testset "Integration with ERGM estimation" begin
        Random.seed!(11)
        net = random_net(12, 25; directed=false, seed=5)

        # Custom terms work inside fit_ergm alongside built-ins
        result = fit_ergm(net, [Edges(), ExampleTerm()])
        @test result.converged
        @test length(result.coefficients) == 2
        @test result.model.formula.terms.names == ["edges", "example"]

        # ERGM's StatsAPI accessors work on fits containing user terms
        @test coef(result) === result.coefficients
        @test stderror(result) === result.std_errors
        @test size(vcov(result)) == (2, 2)

        # The bundled example terms are covariate-only and declare it;
        # unknown user terms keep ERGM's conservative fallback (true)
        @test !ERGM.is_dyad_dependent(ExampleTerm())
        @test !ERGM.is_dyad_dependent(WeightedEdges())
        @test ERGM.is_dyad_dependent(SharedNeighborTerm())

        # ERGM's model-construction validation and attribute materialization
        # cover only its built-in attribute terms: a user term referencing a
        # nonexistent attribute passes through unchanged (and unvalidated)
        model = ERGMModel(ERGMFormula([Edges(), InteractionTerm(:no_such_a,
                                                                :no_such_b)]),
                          net)
        @test model.formula.terms[2] isa InteractionTerm
    end

    @testset "@ergm_term definition checks" begin
        # Complete definition: no warnings
        @test_logs ERGMUserterms._check_term_definition(CompleteTerm)

        # Missing methods produce warnings
        @test_logs (:warn, r"no compute") (:warn, r"no change_stat") (:warn, r"no name") begin
            ERGMUserterms._check_term_definition(IncompleteTerm)
        end
    end

    @testset "Benchmark and profile utilities" begin
        net = random_net(10, 15)
        b = benchmark_term(ExampleTerm(), net; n_iter=10)
        @test b.compute_mean > 0
        @test b.change_stat_mean > 0

        prof = profile_term(ExampleTerm(); sizes=[5, 10], n_iter=5)
        @test length(prof) == 2
        @test prof[1].n_vertices == 5
    end

    @testset "Documentation helpers" begin
        sig = term_signature(TemplateTerm(2.5))
        @test occursin("TemplateTerm", sig)
        @test occursin("param", sig)

        doc = term_documentation(TemplateTerm(2.5))
        @test occursin("template.2.5", doc)
        @test occursin("change_stat", doc)
    end
end
