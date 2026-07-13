using ERGMUserterms
using ERGM
using Networks
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

# ---------------------------------------------------------------------------
# Terms exercising the public term-trait protocol (ERGM.jl src/terms/traits.jl)
# ---------------------------------------------------------------------------

# Declares every trait truthfully: reads a vertex attribute, is defined only on
# directed networks, is dyad-independent, and honours the missing-dyad mask
# (its statistic never reads a masked dyad's face value).
struct HonestTerm <: AbstractUserTerm
    attr::Symbol
end
name(t::HonestTerm) = "honest.$(t.attr)"
function compute(t::HonestTerm, net)
    vals = get_vertex_attribute(net, t.attr)
    total = 0.0
    for e in edges(net)
        i, j = Int(src(e)), Int(dst(e))
        is_missing_dyad(net, i, j) && continue
        total += Float64(get(vals, i, 0.0)) + Float64(get(vals, j, 0.0))
    end
    return total
end
function change_stat(t::HonestTerm, net, i::Int, j::Int)
    is_missing_dyad(net, i, j) && return 0.0
    vals = get_vertex_attribute(net, t.attr)
    return Float64(get(vals, i, 0.0)) + Float64(get(vals, j, 0.0))
end
ERGM.required_vertex_attributes(t::HonestTerm) = (t.attr,)
ERGM.requires_directed(::HonestTerm) = true
ERGM.is_dyad_dependent(::HonestTerm) = false
Networks.supports_missing(::HonestTerm) = true

# Reads a vertex attribute but declares nothing: on a network without :a its
# statistic silently collapses to zero. validate_traits must catch this.
struct UndeclaredAttrTerm <: AbstractUserTerm end
name(::UndeclaredAttrTerm) = "undeclared_attr"
function compute(::UndeclaredAttrTerm, net)
    vals = get_vertex_attribute(net, :a)
    total = 0.0
    for e in edges(net)
        total += Float64(get(vals, Int(src(e)), 0.0))
    end
    return total
end
change_stat(::UndeclaredAttrTerm, net, i::Int, j::Int) =
    Float64(get(get_vertex_attribute(net, :a), i, 0.0))
ERGM.is_dyad_dependent(::UndeclaredAttrTerm) = false

# Counts mutual dyads (genuinely dyad-dependent) but claims independence.
struct LyingDependenceTerm <: AbstractUserTerm end
name(::LyingDependenceTerm) = "lying_dependence"
function compute(::LyingDependenceTerm, net)
    total = 0.0
    for e in edges(net)
        i, j = Int(src(e)), Int(dst(e))
        i < j && has_edge(net, j, i) && (total += 1.0)
    end
    return total
end
change_stat(::LyingDependenceTerm, net, i::Int, j::Int) =
    has_edge(net, j, i) ? 1.0 : 0.0
ERGM.is_dyad_dependent(::LyingDependenceTerm) = false   # false: it is not

# Counts every edge at face value — including masked ones — but claims to
# honour the mask.
struct LyingMissingTerm <: AbstractUserTerm end
name(::LyingMissingTerm) = "lying_missing"
compute(::LyingMissingTerm, net) = Float64(ne(net))
change_stat(::LyingMissingTerm, net, i::Int, j::Int) = 1.0
ERGM.is_dyad_dependent(::LyingMissingTerm) = false
Networks.supports_missing(::LyingMissingTerm) = true    # true: it is not

# The package template shipped in examples/ — included so it is machine-verified
include(joinpath(@__DIR__, "..", "examples", "MyTermPackage", "src",
                 "MyTermPackage.jl"))
using .MyTermPackage: ReciprocatedHomophily

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

        # A user term that declares its attributes (InteractionTerm does, via
        # ERGM.required_vertex_attributes) is validated at model construction
        # exactly like a built-in one
        @test_throws ArgumentError ERGMModel(
            ERGMFormula([Edges(), InteractionTerm(:no_such_a, :no_such_b)]), net)

        # A term declaring no attributes still passes through unchanged
        model = ERGMModel(ERGMFormula([Edges(), ExampleTerm()]), net)
        @test model.formula.terms[2] isa ExampleTerm
    end

    @testset "Public term-trait protocol" begin
        Random.seed!(23)
        net = random_net(10, 22; directed=true, seed=9)
        set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:10))

        # Declarations reach ERGM.jl's public trait functions
        term = HonestTerm(:a)
        @test ERGM.required_vertex_attributes(term) == (:a,)
        @test ERGM.required_edge_attributes(term) == ()
        @test ERGM.requires_directed(term)
        @test !ERGM.requires_undirected(term)
        @test !ERGM.is_dyad_dependent(term)
        @test supports_missing(term)

        # Defaults for a term that declares nothing
        @test ERGM.required_vertex_attributes(ExampleTerm()) == ()
        @test ERGM.required_edge_attributes(ExampleTerm()) == ()
        @test !ERGM.requires_directed(ExampleTerm())
        @test !supports_missing(ExampleTerm())

        # A term declaring all four traits is accepted by ERGMModel and fits
        model = ERGMModel(ERGMFormula([Edges(), term]), net)
        @test model.formula.terms.names == ["edges", "honest.a"]
        @test length(coef(fit_ergm(net, [Edges(), term]))) == 2

        # A declared direction requirement is enforced: rejected undirected
        und = network(10; directed=false)
        set_vertex_attribute!(und, :a, Dict(v => Float64(v) for v in 1:10))
        add_edge!(und, 1, 2)
        @test_throws ArgumentError ERGMModel(ERGMFormula([Edges(), term]), und)

        # A declared vertex attribute the network lacks: the standard error
        @test_throws ArgumentError ERGMModel(
            ERGMFormula([Edges(), HonestTerm(:no_such)]), net)
        err = try
            ERGMModel(ERGMFormula([HonestTerm(:no_such)]), net)
        catch e
            e
        end
        @test occursin("vertex attribute :no_such", err.msg)

        # Backward compatibility: the private names TERGM.jl declares methods
        # on still work, and dispatch on the same generic as the public ones
        @test ERGM._requires_directed === ERGM.requires_directed
        @test ERGM._requires_directed(term)
        @test ERGM._vertex_attribute(term) == :a
        @test ERGM._vertex_attribute(ExampleTerm()) === nothing
    end

    @testset "validate_traits exercises the declarations" begin
        Random.seed!(29)
        net = random_net(10, 24; directed=true, seed=13)
        set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:10))
        set_vertex_attribute!(net, :b, Dict(v => Float64(11 - v) for v in 1:10))

        # Truthful declarations pass (interface + traits)
        @test validate_traits(HonestTerm(:a), net; verbose=false)
        @test validate_term(HonestTerm(:a), net; verbose=false)

        # Reads :a but declares nothing -> caught
        @test !validate_traits(UndeclaredAttrTerm(), net; verbose=false)
        @test !validate_term(UndeclaredAttrTerm(), net; verbose=false)

        # Claims dyad-independence but reads the reverse arc -> caught
        @test !validate_traits(LyingDependenceTerm(), net; verbose=false)

        # Claims to honour the mask but counts masked edges at face value
        @test !validate_traits(LyingMissingTerm(), net; verbose=false)

        # A directed-only term validated on an undirected network is a
        # mismatch the harness reports rather than silently accepting
        und = random_net(10, 20; directed=false, seed=17)
        set_vertex_attribute!(und, :a, Dict(v => Float64(v) for v in 1:10))
        @test !validate_traits(HonestTerm(:a), und; verbose=false)

        # Terms declaring nothing keep validating as before
        @test validate_traits(ExampleTerm(), net; verbose=false)
        @test validate_term(TemplateTerm(1.5), net; verbose=false)

        # Trait checks can be switched off
        @test validate_term(UndeclaredAttrTerm(), net; verbose=false, traits=false)
    end

    @testset "Package template (examples/MyTermPackage)" begin
        net = network(8; directed=true)
        set_vertex_attribute!(net, :group,
                              Dict(v => isodd(v) ? "a" : "b" for v in 1:8))
        for (i, j) in ((1, 3), (3, 1), (2, 4), (4, 2), (1, 2), (5, 7), (6, 8))
            add_edge!(net, i, j)
        end
        term = ReciprocatedHomophily(:group)

        # Statistic and the four declarations
        @test compute(term, net) == 2.0
        @test ERGM.required_vertex_attributes(term) == (:group,)
        @test ERGM.requires_directed(term)
        @test ERGM.is_dyad_dependent(term)
        @test supports_missing(term)

        # The template passes the full harness and ERGM model construction
        Random.seed!(31)
        @test validate_term(term, net; verbose=false)
        @test consistency_check(term, net; exhaustive=true)
        model = ERGMModel(ERGMFormula([Edges(), term]), net)
        @test model.formula.terms.names == ["edges", "recip_homophily.group"]

        # Masked dyads do not enter the statistic (supports_missing = true)
        masked = deepcopy(net)
        set_missing_dyad!(masked, 1, 3)
        @test compute(term, masked) == 1.0
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
