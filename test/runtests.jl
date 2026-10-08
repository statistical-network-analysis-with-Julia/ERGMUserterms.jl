using ERGMUserterms
using ERGM
using NetworkCore
using Graphs
using Logging
using Random
using REPL      # REPL.softscope: the docs-execution gate runs blocks as the REPL would
using TOML
using Test
using Aqua

import ERGMUserterms: name, compute, change_stat

# Text files read by the tests are compared line by line; a Windows checkout
# (git's core.autocrlf) gives them CRLF endings, so normalise to LF.
_readtext(path) = replace(read(path, String), "\r\n" => "\n")

# The workflow testset runs the CI layout step against the sibling checkouts
# that [sources] names. A lone checkout or a registry install (where [sources]
# is not used) has none of them beside it: the step is then not run, with a
# message. Inside the layout every sibling must be present, so a partial
# layout runs the step and fails on the missing one rather than skipping.
function layout_siblings(pkgdir::AbstractString, expected, pkg::AbstractString)
    siblings = sort!([s for s in expected if s != "$pkg.jl"])
    present = filter(s -> isfile(joinpath(dirname(pkgdir), s, "Project.toml")), siblings)
    return (siblings=siblings, present=present, in_layout=!isempty(present))
end

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

# Records every dyad its change_stat is asked about, so the harness's dyad
# sequence can be compared between calls (the rng contract) and counted
struct RecordingTerm <: AbstractUserTerm
    visited::Vector{Tuple{Int, Int}}
end
RecordingTerm() = RecordingTerm(Tuple{Int, Int}[])
name(::RecordingTerm) = "recording"
compute(::RecordingTerm, net) = Float64(ne(net))
function change_stat(t::RecordingTerm, net, i::Int, j::Int)
    push!(t.visited, (i, j))
    return 1.0
end
ERGM.is_dyad_dependent(::RecordingTerm) = false

# Silence the harness's @info/@warn chatter and its println report
quietly(f) = with_logger(f, NullLogger())
silently(f) = redirect_stdout(() -> quietly(f), devnull)

# Extract the `rng=Xoshiro(0x…, …)` replay literal from a warning message and
# turn it back into an rng
function replay_rng(msg::AbstractString)
    m = match(r"rng=(Xoshiro\((?:0x[0-9a-f]{16}(?:, )?){4}\))", msg)
    m === nothing && error("no replay literal in: $msg")
    return eval(Meta.parse(m.captures[1]))
end

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
# The `get(vals, v, 0.0)` zero-fill below is unreachable through `ERGMModel`
# for the DECLARED attribute: `ERGM.Extension.validate_formula` refuses a network on
# which `:attr` is absent or set on only some vertices (statnet refuses NA)
# before any statistic is computed. It exists only for raw `compute` calls.
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
NetworkCore.supports_missing(::HonestTerm) = true

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
NetworkCore.supports_missing(::LyingMissingTerm) = true    # true: it is not

# Defined on undirected networks only (ERGM.requires_undirected), like ERGM's
# own Degree/Kstar/GWDegree
struct UndirectedOnlyTerm <: AbstractUserTerm end
name(::UndirectedOnlyTerm) = "undirected_only"
compute(::UndirectedOnlyTerm, net) = Float64(ne(net))
change_stat(::UndirectedOnlyTerm, net, i::Int, j::Int) = 1.0
ERGM.requires_undirected(::UndirectedOnlyTerm) = true
ERGM.is_dyad_dependent(::UndirectedOnlyTerm) = false

# ... and one whose change_stat refuses a directed network outright, the
# way ERGM's Degree/Kstar/GWDegree do — it cannot be timed on one at all
struct StrictUndirectedTerm <: AbstractUserTerm end
name(::StrictUndirectedTerm) = "strict_undirected"
compute(::StrictUndirectedTerm, net) = Float64(ne(net))
function change_stat(::StrictUndirectedTerm, net, i::Int, j::Int)
    is_directed(net) && error("strict_undirected is defined on undirected networks only")
    return 1.0
end
ERGM.requires_undirected(::StrictUndirectedTerm) = true
ERGM.is_dyad_dependent(::StrictUndirectedTerm) = false

# Records what kind of network each `compute` call saw (directedness, size,
# vertex attributes), so `test_term`'s keywords can be shown to reach every
# network it generates
struct ProbeTerm <: AbstractUserTerm
    seen::Vector{Tuple{Bool, Int, Vector{Symbol}}}
end
ProbeTerm() = ProbeTerm(Tuple{Bool, Int, Vector{Symbol}}[])
name(::ProbeTerm) = "probe"
function compute(t::ProbeTerm, net)
    push!(t.seen, (is_directed(net), Int(nv(net)), sort(list_vertex_attributes(net))))
    return Float64(ne(net))
end
change_stat(::ProbeTerm, net, i::Int, j::Int) = 1.0
ERGM.is_dyad_dependent(::ProbeTerm) = false

# --- terms for the harness-power, g(∅) and edge-attribute testsets ----------
# Mutual-like, with a change statistic that forgets the reciprocation
# whenever i > j (wrong on half of all reciprocated dyads)
struct HalfWrongMutual <: AbstractUserTerm end
name(::HalfWrongMutual) = "halfwrongmutual"
compute(::HalfWrongMutual, net) = compute(Mutual(), net)
change_stat(::HalfWrongMutual, net, i::Int, j::Int) =
    (i < j && has_edge(net, j, i)) ? 1.0 : 0.0

# An isolates count written to the OLD guide's advice: wrong on the empty
# network (g(∅) is n, not 0)
_n_isolates(net) = Float64(count(v -> isempty(inneighbors(net, v)) &&
                                         isempty(outneighbors(net, v)), 1:Int(nv(net))))
# add-direction: an endpoint with no tie other than (i,j) itself stops being
# an isolate
function _isolates_change(net, i, j)
    lone(v, w) = all(==(w), outneighbors(net, v)) && all(==(w), inneighbors(net, v))
    return -Float64(lone(i, j)) - Float64(lone(j, i))
end
struct GuideIsolates <: AbstractUserTerm end
name(::GuideIsolates) = "guideisolates"
compute(::GuideIsolates, net) = (ne(net) == 0 && return 0.0; _n_isolates(net))
change_stat(::GuideIsolates, net, i::Int, j::Int) = _isolates_change(net, i, j)

# ... and without the special case: g(∅) = n, correctly
struct PlainIsolates <: AbstractUserTerm end
name(::PlainIsolates) = "plainisolates"
compute(::PlainIsolates, net) = _n_isolates(net)
change_stat(::PlainIsolates, net, i::Int, j::Int) = _isolates_change(net, i, j)

# The number of NON-ties: dyad-independent, with g(∅) = the number of dyads
# (every built-in dyad-independent term has g(∅) = 0)
struct NonTies <: AbstractUserTerm end
name(::NonTies) = "nonties"
compute(::NonTies, net) = (n = Int(nv(net)); Float64((is_directed(net) ? n * (n - 1) : n * (n - 1) ÷ 2) - ne(net)))
change_stat(::NonTies, net, i::Int, j::Int) = -1.0
ERGM.is_dyad_dependent(::NonTies) = false

# Edge count whose change statistic is wrong on the empty network only
struct WrongWhenEmpty <: AbstractUserTerm end
name(::WrongWhenEmpty) = "wrongwhenempty"
compute(::WrongWhenEmpty, net) = Float64(ne(net))
change_stat(::WrongWhenEmpty, net, i::Int, j::Int) = ne(net) == 0 ? 2.0 : 1.0

# Declares dyad-independence, but the change statistic of dyad (1,2) reads
# two unrelated arcs, (17,18) and (18,17) — whether they DIFFER, so the lie
# is invisible on the empty and on the complete network and only toggling
# one of them at a time exposes it
struct ReadsOneDyad <: AbstractUserTerm end
name(::ReadsOneDyad) = "readsonedyad"
_odd(net) = has_edge(net, 17, 18) != has_edge(net, 18, 17)
compute(::ReadsOneDyad, net) = Float64(ne(net)) + (has_edge(net, 1, 2) && _odd(net) ? 1.0 : 0.0)
change_stat(::ReadsOneDyad, net, i::Int, j::Int) = (i, j) == (1, 2) && _odd(net) ? 2.0 : 1.0
ERGM.is_dyad_dependent(::ReadsOneDyad) = false

# Reads an edge attribute LIVE (no snapshot): fine for every numeric check,
# wrong under MCMC
struct LiveWeight <: AbstractUserTerm end
name(::LiveWeight) = "liveweight"
_live_weight(net, i, j) = (w = get_edge_attribute(net, :weight, i, j);
                           w === nothing ? 1.0 : Float64(w))
function compute(::LiveWeight, net)
    total = 0.0
    for e in edges(net)
        total += _live_weight(net, Int(src(e)), Int(dst(e)))
    end
    return total
end
change_stat(::LiveWeight, net, i::Int, j::Int) = _live_weight(net, i, j)
ERGM.is_dyad_dependent(::LiveWeight) = false

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
                @test consistency_check(term, net; rng=Xoshiro(7))
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
        net = random_net(12, 25; directed=false, seed=5)

        # Custom terms work inside fit_ergm alongside built-ins
        result = fit_ergm(net, [Edges(), ExampleTerm()]; rng=Xoshiro(11))
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
        @test length(coef(fit_ergm(net, [Edges(), term]; rng=Xoshiro(23)))) == 2

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

        # The public trait names are the ones this package and its docs use;
        # ERGM deleted the private `_requires_directed` alias in 0.2.
        @test !isdefined(ERGM, :_requires_directed)
        @test ERGM.requires_directed(term)
        @test ERGM.required_vertex_attributes(term) == (:a,)
        @test ERGM.required_vertex_attributes(ExampleTerm()) == ()
        @test ERGM.has_dyad_dependent(ERGMModel(ERGMFormula([Edges(), term]), net)) == false
    end

    @testset "Declared attributes must be complete (ERGM's NA-completeness rule)" begin
        # `:a` on vertices 1:9 of 10: statnet refuses an NA attribute value,
        # and so does ERGM.Extension.validate_formula — through the public trait
        # protocol, so a user term gets exactly the built-in behaviour
        net = random_net(10, 20; directed=true, seed=41)
        set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:9))
        err = try
            ERGMModel(ERGMFormula([HonestTerm(:a)]), net)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("on every vertex", err.msg)
        @test occursin("vertices 10", err.msg)
        @test_throws ArgumentError ERGMModel(ERGMFormula([HonestTerm(:a)]), net)

        # The harness reports it, in the attribute step and via ERGMModel
        @test !validate_traits(HonestTerm(:a), net; verbose=false, rng=Xoshiro(1))
        @test !validate_term(HonestTerm(:a), net; verbose=false, rng=Xoshiro(1))
        @test_logs (:warn, r"every vertex") match_mode=:any begin
            validate_traits(HonestTerm(:a), net; rng=Xoshiro(1))
        end
        # ... and prints ERGM's own message, not an `ArgumentError(...)` wrapper
        @test_logs (:warn, r"rejected the term or the network: ArgumentError: term") match_mode=:any begin
            validate_traits(HonestTerm(:a), net; rng=Xoshiro(1))
        end

        # Completing the attribute makes the same term valid
        set_vertex_attribute!(net, :a, 10, 10.0)
        @test validate_term(HonestTerm(:a), net; verbose=false, rng=Xoshiro(1))
        @test ERGMModel(ERGMFormula([HonestTerm(:a)]), net) isa ERGMModel
    end

    @testset "validate_traits exercises the declarations" begin
        net = random_net(10, 24; directed=true, seed=13)
        set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:10))
        set_vertex_attribute!(net, :b, Dict(v => Float64(11 - v) for v in 1:10))
        rng() = Xoshiro(29)

        # Truthful declarations pass (interface + traits)
        @test validate_traits(HonestTerm(:a), net; verbose=false, rng=rng())
        @test validate_term(HonestTerm(:a), net; verbose=false, rng=rng())

        # Reads :a but declares nothing -> caught
        @test !validate_traits(UndeclaredAttrTerm(), net; verbose=false, rng=rng())
        @test !validate_term(UndeclaredAttrTerm(), net; verbose=false, rng=rng())

        # Claims dyad-independence but reads the reverse arc -> caught
        @test !validate_traits(LyingDependenceTerm(), net; verbose=false, rng=rng())

        # Claims to honour the mask but counts masked edges at face value
        @test !validate_traits(LyingMissingTerm(), net; verbose=false, rng=rng())

        # A directed-only term validated on an undirected network is a
        # mismatch the harness reports rather than silently accepting
        und = random_net(10, 20; directed=false, seed=17)
        set_vertex_attribute!(und, :a, Dict(v => Float64(v) for v in 1:10))
        @test !validate_traits(HonestTerm(:a), und; verbose=false, rng=rng())

        # Terms declaring nothing keep validating as before
        @test validate_traits(ExampleTerm(), net; verbose=false, rng=rng())
        @test validate_term(TemplateTerm(1.5), net; verbose=false, rng=rng())

        # Trait checks can be switched off
        @test validate_term(UndeclaredAttrTerm(), net; verbose=false, traits=false,
                            rng=rng())

        # A term that does NOT declare supports_missing is validated at the
        # face value of any masked dyad — and the harness says so (the
        # missing-data contract: no face-value number without a word). The
        # verdict is unchanged; the notice is one @info per entry point.
        masked = deepcopy(net)
        set_missing_dyad!(masked, 1, 2)
        notice = r"1 masked dyad; the term declares supports_missing = false, so they are validated at their face value"
        @test_logs (:info, notice) match_mode=:any validate_term(ExampleTerm(), masked; rng=rng())
        @test_logs (:info, notice) match_mode=:any change_stat_check(ExampleTerm(), masked; rng=rng())
        @test_logs (:info, notice) match_mode=:any consistency_check(ExampleTerm(), masked; verbose=true, rng=rng())
        @test validate_term(ExampleTerm(), masked; verbose=false, rng=rng())
        mentions_mask(f) = any(occursin("validated at their face value", l.message)
                               for l in Test.collect_test_logs(f)[1])
        # ... not for a term that honours the mask, and not without masked dyads
        @test !mentions_mask(() -> validate_term(HonestTerm(:a), masked; rng=rng()))
        @test !mentions_mask(() -> validate_term(ExampleTerm(), net; rng=rng()))
        # ... and never with verbose=false (the flag is the whole output)
        @test !mentions_mask(() -> validate_term(ExampleTerm(), masked; verbose=false, rng=rng()))
        @test !mentions_mask(() -> change_stat_check(ExampleTerm(), masked; verbose=false, rng=rng()))
        @test !mentions_mask(() -> consistency_check(ExampleTerm(), masked; rng=rng()))
        # The face value IS what compute reports (8.0 = the masked (1,2) counted)
        @test compute(ExampleTerm(), masked) == compute(ExampleTerm(), net)
    end

    @testset "Two-mode networks are refused, not misjudged" begin
        # On `network(6; bipartite=3)` a within-mode pair cannot hold an edge
        # (`add_edge!` returns false), so the harness's brute-force reference
        # would read 0 there and a CORRECT term would be reported
        # inconsistent (it used to report `Inconsistency at dyad (1,3):
        # change_stat=4.0, brute-force=0.0`). Every entry point that takes a
        # network now refuses a two-mode one up front, with the one-mode-only
        # message ERGMModel itself gives, before any dyad is visited.
        bp = network(6; directed=false, bipartite=3)
        add_edge!(bp, 1, 4); add_edge!(bp, 2, 5)
        @test !add_edge!(bp, 1, 2)
        @test is_two_mode(bp)
        for f in (validate_term, validate_traits, change_stat_check,
                  consistency_check, benchmark_term)
            @test_throws ArgumentError f(ExampleTerm(), bp)
            err = try
                f(ExampleTerm(), bp)
            catch e
                e
            end
            @test startswith(err.msg, "$(nameof(f)) validates terms on networks ERGMModel accepts only")
            # ... and the reason is ERGM's own sentence
            # (`ERGM.Extension.require_supported_network`), not a copy kept in
            # step by hand
            @test occursin("this network is two-mode (bipartite)", err.msg)
            @test occursin("ERGM.jl fits one-mode networks only", err.msg)
            ergm_err = try ERGM.Extension.require_supported_network(bp) catch e; e end
            @test endswith(err.msg, ergm_err.msg)
        end
        @test_throws ArgumentError consistency_check(ExampleTerm(), bp; exhaustive=true)
        # ... which is exactly the model's own verdict
        @test_throws ArgumentError ERGMModel(ERGMFormula([Edges(), ExampleTerm()]), bp)
        # No dyad is visited and nothing is logged before the refusal
        t = RecordingTerm()
        @test_throws ArgumentError validate_term(t, bp)
        @test isempty(t.visited)
        @test_logs (try validate_term(ExampleTerm(), bp) catch; end)
        # A one-mode network of the same size is validated as before, and
        # test_term/profile_term build one-mode networks of their own
        @test validate_term(ExampleTerm(), network(6; directed=false); verbose=false, rng=Xoshiro(1))
        @test silently(() -> test_term(ExampleTerm(); n_vertices=6, n_tests=3, rng=Xoshiro(1)))
    end

    @testset "Self-loop networks are refused up front, not blamed on the term" begin
        # ERGM.jl refuses a network containing a self-loop (its statistics
        # would count the loop; its estimators and samplers never touch the
        # diagonal). The harness used to validate such a network dyad by dyad — every
        # check green, `compute` counting the loop — and only failed at the
        # ERGMModel step with a message that read as a defect of the term,
        # while `traits=false`, `change_stat_check` and `consistency_check`
        # returned `true` outright. Now every entry point refuses it first,
        # with ERGM's own reasoning, and no dyad is visited.
        nl = network(6; directed=true, loops=true)
        add_edge!(nl, 1, 1); add_edge!(nl, 1, 2); add_edge!(nl, 2, 3)
        @test has_edge(nl, 1, 1)
        @test_throws ArgumentError ERGMModel(ERGMFormula([Edges(), ExampleTerm()]), nl)
        for f in (validate_term, validate_traits, change_stat_check,
                  consistency_check, benchmark_term)
            @test_throws ArgumentError f(ExampleTerm(), nl)
            err = try
                f(ExampleTerm(), nl)
            catch e
                e
            end
            @test startswith(err.msg, "$(nameof(f)) validates terms on networks ERGMModel accepts only")
            @test occursin("contains 1 self-loop (at vertex 1)", err.msg)
            @test occursin("ERGM.jl models the off-diagonal dyads only", err.msg)
            @test occursin("rem_edge!(net, v, v)", err.msg)
            # ... ERGM's own sentence (`ERGM.Extension.require_supported_network`)
            ergm_err = try ERGM.Extension.require_supported_network(nl) catch e; e end
            @test endswith(err.msg, ergm_err.msg)
        end
        @test_throws ArgumentError validate_term(ExampleTerm(), nl; traits=false)
        @test_throws ArgumentError consistency_check(ExampleTerm(), nl; exhaustive=true)
        t = RecordingTerm()
        @test_throws ArgumentError validate_term(t, nl)
        @test isempty(t.visited)
        @test_logs (try validate_term(ExampleTerm(), nl) catch; end)
        # Several loops: all listed
        add_edge!(nl, 4, 4)
        err = try; change_stat_check(ExampleTerm(), nl); catch e; e; end
        @test occursin("2 self-loops (at vertex 1, 4)", err.msg)
        # A loops-ALLOWED network without any loop present is fine
        # (`loops=true` is a capability; ERGM refuses a loop, not the flag)
        rem_edge!(nl, 1, 1); rem_edge!(nl, 4, 4)
        @test validate_term(ExampleTerm(), nl; verbose=false, rng=Xoshiro(1))
        @test consistency_check(ExampleTerm(), nl; exhaustive=true)
    end

    @testset "No verdict without evidence" begin
        # A PASS resting on zero checked dyads is the defect the harness
        # exists to catch. `n_tests=0` ("skip the only numeric check"), a
        # network with fewer than 2 vertices (no dyad exists) and an
        # out-of-range `density` are ArgumentErrors on every entry point,
        # consistent with what `benchmark_term` always did.
        net = random_net(10, 15; seed=3)
        for f in (validate_term, change_stat_check)
            @test_throws ArgumentError f(ToggleConventionTerm(), net; n_tests=0)
            @test_throws ArgumentError f(ExampleTerm(), net; n_tests=0)
        end
        for f in (validate_term, validate_traits, change_stat_check,
                  consistency_check, benchmark_term), n in (0, 1)
            @test_throws ArgumentError f(ExampleTerm(), network(n; directed=true))
            err = try; f(ExampleTerm(), network(n; directed=true)); catch e; e; end
            @test occursin("at least 2 vertices", err.msg)
        end
        @test_throws ArgumentError validate_term(ExampleTerm(), network(1; directed=true); traits=false)
        @test_throws ArgumentError consistency_check(ExampleTerm(), network(1; directed=false); exhaustive=true)
        for kw in ((n_vertices=1,), (n_vertices=0,), (density=1.5,), (density=-0.2,),
                   (n_tests=0,), (n_tests=-3,))
            @test_throws ArgumentError silently(() -> test_term(ExampleTerm(); rng=Xoshiro(1), kw...))
        end
        @test_throws ArgumentError profile_term(ExampleTerm(); sizes=[5], density=1.5, n_iter=2)
        @test_throws ArgumentError profile_term(ExampleTerm(); sizes=[1], n_iter=2)
        @test_throws ArgumentError profile_term(ExampleTerm(); sizes=[5], n_iter=0)
        # the boundaries themselves are fine: density 0 and 1 are proportions
        @test silently(() -> test_term(ExampleTerm(); n_vertices=4, density=0.0, n_tests=2, rng=Xoshiro(1)))
        @test silently(() -> test_term(ExampleTerm(); n_vertices=4, density=1.0, n_tests=2, rng=Xoshiro(1)))
        @test silently(() -> test_term(ExampleTerm(); n_vertices=2, n_tests=1, rng=Xoshiro(1)))
        @test silently(() -> test_term(ExampleTerm(); n_vertices=2, n_tests=1, directed=false, rng=Xoshiro(1)))

        # A toggle-direction term cannot pass through ANY accepted keyword
        # combination of any validator: the smallest accepted run still
        # checks one dyad in both states
        for n_tests in (1, 2, 5), traits in (true, false), n in (2, 3, 6), directed in (true, false)
            g = random_net(n, directed ? n : 1; directed=directed, seed=n)
            @test !validate_term(ToggleConventionTerm(), g; verbose=false, traits=traits,
                                 n_tests=n_tests, rng=Xoshiro(n_tests))
            @test !change_stat_check(ToggleConventionTerm(), g; verbose=false,
                                     n_tests=n_tests, rng=Xoshiro(n_tests))
        end
        for n in (2, 3, 6), directed in (true, false), exhaustive in (true, false)
            g = random_net(n, directed ? n : 1; directed=directed, seed=n)
            @test !consistency_check(ToggleConventionTerm(), g; exhaustive=exhaustive, rng=Xoshiro(n))
        end
        for n in (2, 3, 8), directed in (true, false), density in (0.0, 0.3, 1.0)
            @test !silently(() -> test_term(ToggleConventionTerm(); n_vertices=n, density=density,
                                            n_tests=1, directed=directed, rng=Xoshiro(n)))
        end
    end

    @testset "A constant attribute is reported as undecided, never as ✓" begin
        # Perturbing a constant attribute changes nothing, so whether the
        # term reads it cannot be decided; the harness used to skip the check and
        # still printed "✓ declared attributes match the ones the term
        # reads" — letting an undeclared-attribute term (the silent all-zero
        # column trap) pass with a positive finding it never made.
        cn = random_net(8, 12; directed=true, seed=5)
        set_vertex_attribute!(cn, :a, Dict(v => 1.0 for v in 1:8))
        undecidable = r"attribute :a is constant on this network, so whether the term reads it is undecidable; validate on a network where it varies"
        logs = Test.collect_test_logs(() -> validate_traits(UndeclaredAttrTerm(), cn; rng=Xoshiro(1)))[1]
        @test any(l.level == Logging.Warn && occursin(undecidable, l.message) for l in logs)
        @test !any(occursin("✓ declared attributes match", l.message) for l in logs)
        @test any(occursin("except that :a could not be tested (constant on this network)", l.message)
                  for l in logs)
        # Undeclared: no evidence of a violation, so the verdict stands
        # (with the warning); declared: the network is not one the term is
        # meant for, so the run FAILS
        @test validate_traits(UndeclaredAttrTerm(), cn; verbose=false, rng=Xoshiro(1))
        @test !validate_traits(HonestTerm(:a), cn; verbose=false, rng=Xoshiro(1))
        @test !validate_term(HonestTerm(:a), cn; verbose=false, rng=Xoshiro(1))
        @test_logs (:warn, r"declares required vertex attribute :a, which is constant on this network") match_mode=:any begin
            validate_traits(HonestTerm(:a), cn; rng=Xoshiro(1))
        end
        # Make it vary and the same terms get their real verdicts
        set_vertex_attribute!(cn, :a, Dict(v => Float64(v) for v in 1:8))
        @test !validate_traits(UndeclaredAttrTerm(), cn; verbose=false, rng=Xoshiro(1))
        @test validate_traits(HonestTerm(:a), cn; verbose=false, rng=Xoshiro(1))
        @test_logs (:info, r"✓ declared attributes match the ones the term reads") match_mode=:any begin
            validate_traits(HonestTerm(:a), cn; rng=Xoshiro(1))
        end
        # A constant EDGE attribute gets the same treatment
        add_edge!(cn, 1, 2)
        for e in edges(cn)
            set_edge_attribute!(cn, :w, src(e), dst(e), 2.0)
        end
        @test_logs (:warn, r"edge attribute :w is constant on this network") match_mode=:any begin
            validate_traits(WeightedEdges(:w), cn; rng=Xoshiro(1))
        end
    end

    @testset "consistency_check visits each unordered dyad once on an undirected network" begin
        # It used to enumerate ordered pairs regardless of directedness: on a
        # 5-vertex undirected network the 10 dyads cost 40 change_stat calls
        # (each dyad twice, two calls per check). Now (i, j) with j < i is
        # skipped, as `_term_fingerprint` does, and the random budget is
        # sized by `ERGM.Extension.n_observed_dyads`.
        u = network(5; directed=false)
        add_edge!(u, 1, 2); add_edge!(u, 3, 4)
        t = RecordingTerm()
        @test consistency_check(t, u; exhaustive=true)
        @test length(t.visited) == 2 * 10          # as-is + toggled, per dyad
        seen = t.visited[1:2:end]
        @test t.visited[1:2:end] == t.visited[2:2:end]
        @test seen == [(i, j) for i in 1:5 for j in 1:5 if i < j]
        # random mode: canonical (min, max) dyads, no repeats, budget = n(n-1)/2
        t = RecordingTerm()
        @test consistency_check(t, u; rng=Xoshiro(4))
        seen = t.visited[1:2:end]
        @test all(i < j for (i, j) in seen)
        @test allunique(seen)
        @test length(seen) <= 10
        # directed: every ordered pair, as before
        d = network(5; directed=true)
        add_edge!(d, 1, 2)
        t = RecordingTerm()
        @test consistency_check(t, d; exhaustive=true)
        @test length(t.visited) == 2 * 20
        @test allunique(t.visited[1:2:end])
        # the verdict is unaffected on either kind
        @test !consistency_check(ToggleConventionTerm(), u; exhaustive=true)
        @test !consistency_check(ToggleConventionTerm(), u; rng=Xoshiro(4))
    end

    @testset "profile_term builds the network kinds test_term can" begin
        # An undirected-only term whose change_stat errors on a directed
        # network could not be profiled at all (the profiler built directed
        # networks unconditionally), and an attribute-declaring term was
        # timed on its attribute-absent fallback branch
        @test_throws ErrorException profile_term(StrictUndirectedTerm(); sizes=[5], n_iter=2, rng=Xoshiro(1))
        @test length(profile_term(StrictUndirectedTerm(); sizes=[5], n_iter=2, directed=false,
                                  rng=Xoshiro(1))) == 1
        probe = ProbeTerm()
        prof = profile_term(probe; sizes=[5, 8], n_iter=2, directed=false,
                            vertex_attributes=Dict(:x => (v -> Float64(v))), rng=Xoshiro(1))
        @test [r.n_vertices for r in prof] == [5, 8]
        @test !isempty(probe.seen)
        @test all(!d for (d, _, _) in probe.seen)
        @test all(attrs == [:x] for (_, _, attrs) in probe.seen)
        @test Set(n for (_, n, _) in probe.seen) == Set([5, 8])
        probe = ProbeTerm()
        profile_term(probe; sizes=[5], n_iter=2, rng=Xoshiro(1))
        @test all(d for (d, _, _) in probe.seen)
        @test all(isempty(attrs) for (_, _, attrs) in probe.seen)
        # an undirected network of n vertices has at most n(n-1)/2 edges
        prof = profile_term(TemplateTerm(1.0); sizes=[6], density=1.0, n_iter=2,
                            directed=false, rng=Xoshiro(1))
        @test prof[1].ne <= 15
        # the attribute really reaches the timed term: InteractionTerm
        # reads :a/:b; a short vector spec is refused
        @test length(profile_term(InteractionTerm(:a, :b); sizes=[4, 6], n_iter=2, rng=Xoshiro(1),
                                  vertex_attributes=Dict(:a => collect(1.0:6.0), :b => collect(6.0:-1:1.0)))) == 2
        @test_throws ArgumentError profile_term(InteractionTerm(:a, :b); sizes=[4, 6], n_iter=2, rng=Xoshiro(1),
                                                vertex_attributes=Dict(:a => [1.0], :b => [1.0]))
        @test issubset([:directed, :vertex_attributes, :sizes, :density, :n_iter, :rng],
                       Base.kwarg_decl(only(methods(profile_term))))
        # equal rngs, equal networks, with the new keywords too
        p1 = profile_term(ExampleTerm(); sizes=[6, 9], n_iter=2, directed=false,
                          vertex_attributes=Dict(:x => (v -> v)), rng=Xoshiro(5))
        p2 = profile_term(ExampleTerm(); sizes=[6, 9], n_iter=2, directed=false,
                          vertex_attributes=Dict(:x => (v -> v)), rng=Xoshiro(5))
        @test [r.ne for r in p1] == [r.ne for r in p2]
    end

    @testset "Bare `using` gets an actionable message" begin
        # The most-documented mistake: methods defined after `using ERGM`
        # (no `import ERGM: …`) are LOCAL functions the shared generics never
        # dispatch to. The harness used to report only ERGM's fallback text
        # (`compute() not implemented for MyTerm`) while the author was
        # looking at the compute method they had just written.
        bare = Module(:BareUsing)
        Core.eval(bare, quote
            using ERGM, ERGMUserterms, NetworkCore
            struct MyTerm <: AbstractUserTerm end
            name(::MyTerm) = "myterm"
            compute(::MyTerm, net) = Float64(ne(net))
            change_stat(::MyTerm, net, i::Int, j::Int) = 1.0
        end)
        bare_term = Core.eval(bare, :(MyTerm()))
        @test Core.eval(bare, :(compute)) !== ERGM.compute        # the trap, reproduced
        net = random_net(6, 8; seed=1)
        import_hint = r"if you defined one after a bare `using`, it is a local function ERGM.jl never calls: add `import ERGM: name, compute, change_stat` before the definitions"
        @test !validate_term(bare_term, net; verbose=false, rng=Xoshiro(1))
        @test_logs (:info, r"name\(\) returns: .*MyTerm \(ERGM's fallback — no name method for MyTerm reaches the shared generic") (:warn, r"compute\(\) failed: ArgumentError: compute\(\) not implemented for .*MyTerm.* — no compute method for MyTerm reaches the shared generic") (:info, r"Testing change_stat") (:warn, r"change_stat\(\d+, \d+\) failed: ArgumentError: change_stat\(\) not implemented for .*MyTerm.* — no change_stat method for MyTerm reaches the shared generic.*reproduce with rng=") begin
            validate_term(bare_term, net; rng=Xoshiro(1))
        end
        logs = Test.collect_test_logs(() -> validate_term(bare_term, net; rng=Xoshiro(1)))[1]
        @test count(occursin(import_hint, l.message) for l in logs) == 3
        # @ergm_term's definition check says the same
        @test_logs (:warn, r"no compute\(term, net\) method — no compute method for MyTerm reaches the shared generic.*import ERGM: name, compute, change_stat") (:warn, r"no change_stat.*import ERGM") (:warn, r"no name.*import ERGM") begin
            ERGMUserterms._check_term_definition(typeof(bare_term))
        end
        # A term whose methods DO reach the generics gets no such hint,
        # even when compute fails for another reason
        @test_logs ERGMUserterms._check_term_definition(CompleteTerm)
        logs = Test.collect_test_logs(() -> validate_term(CompleteTerm(), net; rng=Xoshiro(1)))[1]
        @test !any(occursin("bare `using`", l.message) for l in logs)
    end

    # ------------------------------------------------------------------
    # The rng contract (all randomness flows from the caller's rng)
    # ------------------------------------------------------------------
    @testset "rng contract: every entry point is a pure function of rng" begin
        net = random_net(12, 30; directed=true, seed=19)
        set_vertex_attribute!(net, :a, Dict(v => Float64(v) for v in 1:12))

        # Two calls with equal rngs: identical verdict AND identical dyad
        # sequence, for each of the seven entry points
        function visited(f)
            t = RecordingTerm()
            r = silently(() -> f(t))
            return r, copy(t.visited)
        end
        for f in (t -> validate_term(t, net; rng=Xoshiro(5)),
                  t -> validate_traits(t, net; rng=Xoshiro(5)),
                  t -> change_stat_check(t, net; n_tests=6, rng=Xoshiro(5)),
                  t -> consistency_check(t, net; rng=Xoshiro(5)),
                  t -> test_term(t; n_vertices=9, density=0.2, n_tests=5, rng=Xoshiro(5)))
            (r1, v1), (r2, v2) = visited(f), visited(f)
            @test r1 === r2 === true
            @test v1 == v2
            @test !isempty(v1)
            @test all(i != j for (i, j) in v1)          # never a self-dyad
        end
        # benchmark_term/profile_term: timings differ, the dyads may not
        (_, v1), (_, v2) = visited(t -> benchmark_term(t, net; n_iter=20, rng=Xoshiro(5))),
                           visited(t -> benchmark_term(t, net; n_iter=20, rng=Xoshiro(5)))
        @test v1 == v2 && length(v1) == 20
        p1, (_, w1) = let t = RecordingTerm()
            profile_term(t; sizes=[6, 9], n_iter=4, rng=Xoshiro(5)), (nothing, copy(t.visited))
        end
        p2, (_, w2) = let t = RecordingTerm()
            profile_term(t; sizes=[6, 9], n_iter=4, rng=Xoshiro(5)), (nothing, copy(t.visited))
        end
        @test w1 == w2 && length(w1) == 8
        @test [r.ne for r in p1] == [r.ne for r in p2]

        # Different rngs visit different dyads (the keyword is really used)
        (_, a), (_, b) = visited(t -> change_stat_check(t, net; n_tests=6, rng=Xoshiro(5))),
                         visited(t -> change_stat_check(t, net; n_tests=6, rng=Xoshiro(6)))
        @test a != b

        # The caller's rng is advanced (draws come from it, not from a copy)
        r = Xoshiro(5)
        change_stat_check(ExampleTerm(), net; n_tests=3, rng=r)
        @test r != Xoshiro(5)

        # The keywords are declared on every entry point
        for f in (validate_term, validate_traits, change_stat_check,
                  consistency_check, test_term, benchmark_term, profile_term)
            @test :rng in Base.kwarg_decl(only(methods(f)))
        end
        @test :n_tests in Base.kwarg_decl(only(methods(validate_term)))
        @test :n_tests in Base.kwarg_decl(only(methods(test_term)))
        @test issubset([:directed, :vertex_attributes],
                       Base.kwarg_decl(only(methods(test_term))))
    end

    @testset "Failures carry a replayable rng literal" begin
        net = random_net(10, 18; directed=true, seed=21)
        lit = r"rng=Xoshiro\(0x[0-9a-f]{16}(, 0x[0-9a-f]{16}){3}\)"

        # change_stat_check: the state-independence failure names the dyad
        # and the rng; replaying the literal visits the very same dyad
        @test_logs (:warn, lit) match_mode=:any begin
            change_stat_check(ToggleConventionTerm(), net; n_tests=1, rng=Xoshiro(3))
        end
        msg = Test.collect_test_logs(() ->
            change_stat_check(ToggleConventionTerm(), net; n_tests=1, rng=Xoshiro(3)))[1][1].message
        @test occursin("passed to change_stat_check", msg)
        dyad = match(r"dyad \((\d+),(\d+)\)", msg)
        t = RecordingTerm()
        quietly(() -> change_stat_check(t, net; n_tests=1, rng=replay_rng(msg)))
        @test t.visited[1] == (parse(Int, dyad[1]), parse(Int, dyad[2]))

        # validate_traits: the supports_missing failure carries it too
        @test_logs (:warn, lit) match_mode=:any begin
            validate_traits(LyingMissingTerm(), net; rng=Xoshiro(3))
        end
        msg = Test.collect_test_logs(() ->
            validate_traits(LyingMissingTerm(), net; rng=Xoshiro(3)))[1]
        w = only(l for l in msg if l.level == Logging.Warn)
        @test occursin("passed to validate_traits", w.message)
        # ... and the dyad-dependence failure
        w = only(l for l in Test.collect_test_logs(() ->
                     validate_traits(LyingDependenceTerm(), net; rng=Xoshiro(3)))[1]
                 if l.level == Logging.Warn)
        @test occursin(lit, w.message)

        # validate_term's own failure line and consistency_check's
        @test_logs (:warn, lit) match_mode=:any begin
            validate_term(ToggleConventionTerm(), net; rng=Xoshiro(3))
        end
        @test_logs (:warn, lit) match_mode=:any begin
            consistency_check(ToggleConventionTerm(), net; verbose=true, rng=Xoshiro(3))
        end
        # consistency_check is silent unless asked
        @test_logs consistency_check(ToggleConventionTerm(), net; rng=Xoshiro(3))

        # The default rng (TaskLocalRNG) gets a literal too; a MersenneTwister
        # has no short literal, so the message simply has no hint
        w = only(l for l in Test.collect_test_logs(() ->
                     change_stat_check(ToggleConventionTerm(), net; n_tests=1))[1]
                 if l.level == Logging.Warn)
        @test occursin(lit, w.message)
        w = only(l for l in Test.collect_test_logs(() ->
                     change_stat_check(ToggleConventionTerm(), net; n_tests=1,
                                       rng=MersenneTwister(1)))[1]
                 if l.level == Logging.Warn)
        @test !occursin("reproduce", w.message)

        # test_term prints the literal up front so a failure inside its
        # generated network is replayable from the report alone
        out = mktemp() do path, io
            redirect_stdout(io) do
                quietly(() -> test_term(ToggleConventionTerm(); n_vertices=8, rng=Xoshiro(2)))
            end
            flush(io)
            _readtext(path)
        end
        @test occursin(lit, out)
        @test occursin("Some tests FAILED", out)
    end

    @testset "test_term runs and honours n_tests" begin
        @test silently(() -> test_term(ExampleTerm(); n_vertices=8, density=0.2,
                                       n_tests=4, rng=Xoshiro(2)))
        @test !silently(() -> test_term(ToggleConventionTerm(); rng=Xoshiro(2)))

        # n_tests reaches validate_term: with n_tests=3 the Real-check loop
        # and the consistency check each visit 3 dyads (2 change_stat calls
        # per consistency dyad), and the path check from the empty network
        # one per edge
        t = RecordingTerm()
        silently(() -> test_term(t; n_vertices=8, density=0.2, n_tests=3, rng=Xoshiro(2)))
        # validate_term: 3 (Real check) + 3*2 (consistency) + ne (path) +
        # fingerprints (validate_traits) — pin the interface part by running
        # it alone
        t = RecordingTerm()
        net = random_net(8, 12; seed=2)
        quietly(() -> validate_term(t, net; n_tests=3, traits=false, rng=Xoshiro(2)))
        @test length(t.visited) == 3 + 3 * 2 + ne(net)
        @test allunique(t.visited[1:3])                  # distinct dyads
        # The default is every dyad: 56 + 56*2 + ne
        t = RecordingTerm()
        quietly(() -> validate_term(t, net; traits=false))
        @test length(t.visited) == 56 + 56 * 2 + ne(net)
        @test Set(t.visited[1:56]) == Set((i, j) for i in 1:8 for j in 1:8 if i != j)

        # change_stat_check: exactly n_tests dyads, each seen twice
        # (as-is, then with the dyad toggled)
        t = RecordingTerm()
        @test change_stat_check(t, net; n_tests=7, rng=Xoshiro(2))
        @test length(t.visited) == 14
        @test t.visited[1:2:end] == t.visited[2:2:end]   # same dyad, both states
        @test allunique(t.visited[1:2:end])              # sampled without replacement
        # ... and every dyad by default
        t = RecordingTerm()
        @test change_stat_check(t, net)
        @test Set(t.visited) == Set((i, j) for i in 1:8 for j in 1:8 if i != j)
        @test length(t.visited) == 2 * 56

        # n_tests is capped at the number of dyads; 0 (and anything below)
        # is refused — it would mean "skip the only numeric check" (see the
        # "no verdict without evidence" testset)
        tiny = network(3; directed=false); add_edge!(tiny, 1, 2)
        t = RecordingTerm()
        @test quietly(() -> validate_term(t, tiny; n_tests=50, traits=false, rng=Xoshiro(2)))
        @test length(t.visited) == 3 + 3 * 2 + 1
        for bad in (0, -1)
            @test_throws ArgumentError validate_term(ExampleTerm(), tiny; n_tests=bad)
            @test_throws ArgumentError change_stat_check(ExampleTerm(), tiny; n_tests=bad)
        end

        # consistency_check's random mode never draws a self-dyad and never
        # repeats a dyad
        t = RecordingTerm()
        @test consistency_check(t, net; exhaustive=false, rng=Xoshiro(2))
        seen = t.visited[1:2:end]
        @test all(i != j for (i, j) in seen)
        @test allunique(seen)
        @test length(seen) <= min(100, 8 * 7)
    end

    @testset "test_term: directed= and vertex_attributes= reach every generated network" begin
        # The tutorial's step 3 must be able to run the package's own
        # template term, which declares a vertex attribute — and a term that
        # declares requires_undirected. Before 0.2.0 test_term could build
        # only a bare directed network, so both "failed" the tutorial.
        recip = ReciprocatedHomophily(:group)
        @test silently(() -> test_term(recip; n_vertices=8, n_tests=5, rng=Xoshiro(1),
                                       vertex_attributes=Dict(:group => (v -> isodd(v) ? "a" : "b"))))
        # ... and without the attribute the mismatch is still reported, not
        # silently passed (the term declares it; the network lacks it)
        @test !silently(() -> test_term(recip; n_vertices=8, n_tests=5, rng=Xoshiro(1)))
        @test_logs (:warn, r"required vertex attribute :group, which the validation network does not have") match_mode=:any begin
            redirect_stdout(devnull) do
                test_term(recip; n_vertices=8, n_tests=5, rng=Xoshiro(1))
            end
        end

        # A vector spec (one value per vertex) works too; a short one is refused
        inter = InteractionTerm(:a, :b)
        @test silently(() -> test_term(inter; n_vertices=8, n_tests=5, rng=Xoshiro(2),
                                       vertex_attributes=Dict(:a => collect(1.0:8.0),
                                                              :b => collect(8.0:-1:1.0))))
        @test_throws ArgumentError silently(() -> test_term(inter; n_vertices=8,
                                                            vertex_attributes=Dict(:a => [1.0], :b => [2.0])))

        # requires_undirected: testable on an undirected network, and the
        # default directed network is (correctly) a failure
        @test silently(() -> test_term(UndirectedOnlyTerm(); n_vertices=8, n_tests=5,
                                       directed=false, rng=Xoshiro(2)))
        @test !silently(() -> test_term(UndirectedOnlyTerm(); n_vertices=8, n_tests=5,
                                        rng=Xoshiro(2)))

        # Every network test_term builds — random, empty, complete — has the
        # requested directedness and attributes
        probe = ProbeTerm()
        @test silently(() -> test_term(probe; n_vertices=12, n_tests=4, directed=false,
                                       vertex_attributes=Dict(:x => (v -> Float64(v))),
                                       rng=Xoshiro(3)))
        @test !isempty(probe.seen)
        @test all(!d for (d, _, _) in probe.seen)
        @test all(attrs == [:x] for (_, _, attrs) in probe.seen)
        @test Set(n for (_, n, _) in probe.seen) == Set([12, 10])     # random/empty: 12; complete: min(10, 12)
        probe = ProbeTerm()
        silently(() -> test_term(probe; n_vertices=8, n_tests=4, rng=Xoshiro(3)))
        @test all(d for (d, _, _) in probe.seen)
        @test all(isempty(attrs) for (_, _, attrs) in probe.seen)

        # Equal rngs, equal runs, with the new keywords too
        a = let t = RecordingTerm()
            silently(() -> test_term(t; n_vertices=9, n_tests=5, directed=false,
                                     vertex_attributes=Dict(:x => (v -> v)), rng=Xoshiro(5)))
            copy(t.visited)
        end
        b = let t = RecordingTerm()
            silently(() -> test_term(t; n_vertices=9, n_tests=5, directed=false,
                                     vertex_attributes=Dict(:x => (v -> v)), rng=Xoshiro(5)))
            copy(t.visited)
        end
        @test a == b && !isempty(a)
    end

    @testset "No bare random draws in src" begin
        src = _readtext(joinpath(@__DIR__, "..", "src", "ERGMUserterms.jl"))
        @test !occursin(r"\brand\((?!rng\b)", src)
        @test !occursin(r"\brandperm\((?!rng\b)", src)
        @test !occursin(r"\bseed!\(", src)
        @test occursin("rng::AbstractRNG=Random.default_rng()", src)
    end

    @testset "Docstring examples of the harness execute as written" begin
        # Every ```julia block in the seven entry points' docstrings runs in
        # a fresh module (the same rule tools/check_snippets.jl applies to
        # README/docs), so the documented calls cannot rot
        function code_blocks(f)
            multidoc = Base.Docs.meta(ERGMUserterms)[Base.Docs.Binding(ERGMUserterms, nameof(f))]
            text = join((join(String.(d.text), "\n") for d in values(multidoc.docs)), "\n")
            return [m.captures[1] for m in eachmatch(r"```julia\n(.*?)\n```"s, text)]
        end
        for f in (validate_term, validate_traits, change_stat_check,
                  consistency_check, test_term, benchmark_term, profile_term)
            blocks = code_blocks(f)
            @test !isempty(blocks)
            for block in blocks
                m = Module()
                @test (silently(() -> Core.eval(m, Meta.parseall(block))); true)
            end
        end
    end

    @testset "Shared verbs are the shared generics" begin
        @test ERGMUserterms.compute === NetworkCore.compute === ERGM.compute
        @test ERGMUserterms.name === NetworkCore.name === ERGM.name
        @test ERGMUserterms.change_stat === ERGM.change_stat
        @test ERGM.supports_missing === NetworkCore.supports_missing
        @test ERGMUserterms.supports_missing === NetworkCore.supports_missing

        # Fresh process: co-loading the three packages leaves every shared
        # verb defined (no ambiguous-export UndefVarError), including the
        # re-exported interface generics user code extends
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(@__DIR__)) -e 'using ERGM, ERGMUserterms, NetworkCore; compute; name; change_stat; supports_missing; validate_term; @assert compute === NetworkCore.compute; @assert change_stat === ERGM.change_stat; @assert validate_term(ExampleTerm(), network(4; directed=true); verbose=false)'`
        @test success(cmd)
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
        @test validate_term(term, net; verbose=false, rng=Xoshiro(31))
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
        b = benchmark_term(ExampleTerm(), net; n_iter=10, rng=Xoshiro(1))
        @test b.compute_mean > 0
        @test b.change_stat_mean > 0
        @test_throws ArgumentError benchmark_term(ExampleTerm(), network(1); n_iter=2)
        @test_throws ArgumentError benchmark_term(ExampleTerm(), net; n_iter=0)

        prof = profile_term(ExampleTerm(); sizes=[5, 10], n_iter=5, rng=Xoshiro(1))
        @test length(prof) == 2
        @test prof[1].n_vertices == 5
        @test prof[1].ne == profile_term(ExampleTerm(); sizes=[5], n_iter=5, rng=Xoshiro(1))[1].ne
    end

    @testset "Documentation helpers" begin
        sig = term_signature(TemplateTerm(2.5))
        @test occursin("TemplateTerm", sig)
        @test occursin("param", sig)

        doc = term_documentation(TemplateTerm(2.5))
        @test occursin("template.2.5", doc)
        @test occursin("change_stat", doc)
    end

    # ------------------------------------------------------------------
    # Documentation gates
    # ------------------------------------------------------------------
    @testset "Every exported docstring carries a runnable example" begin
        # A docs build with checkdocs=:exports checks presence, not content:
        # walk the docsystem instead. Every export must carry an
        # ERGMUserterms-OWNED docstring (`Base.Docs.meta(ERGMUserterms)` —
        # this includes the docstrings the package attaches to the
        # re-exported generics `name`/`compute`/`change_stat`, whose bindings
        # resolve to ERGM/NetworkCore) with at least one fenced ```julia block,
        # and every such block must RUN in a fresh module that has done
        # nothing but `using ERGMUserterms`, so an example needing
        # `NetworkCore`, `Random` or `Graphs: src, dst` says so itself. Sketches
        # that are deliberately not standalone use a ```jl fence. Mirrors
        # ERGM.jl's testset of the same name.
        meta = Base.Docs.meta(ERGMUserterms)
        documented_elsewhere(b) = any(haskey(Base.Docs.meta(m), b)
                                      for m in (ERGM, NetworkCore, Graphs))
        undocumented = String[]
        missing_example = String[]
        blocks = Tuple{String, String}[]
        for nm in names(ERGMUserterms)
            nm === :ERGMUserterms && continue
            b = Base.Docs.Binding(ERGMUserterms, nm)
            if !haskey(meta, b)
                documented_elsewhere(b) || push!(undocumented, string(nm))
                continue
            end
            has_example = false
            for (_, ds) in meta[b].docs
                txt = ds.text isa AbstractString ? ds.text : join(string.(ds.text), "\n")
                for m in eachmatch(r"```julia\n(.*?)```"s, txt)
                    has_example = true
                    push!(blocks, (string(nm), String(m.captures[1])))
                end
                occursin("```jldoctest", txt) && (has_example = true)
            end
            has_example || push!(missing_example, string(nm))
        end
        @test isempty(undocumented)
        @test isempty(missing_example)
        # The three re-exported generics carry ERGMUserterms-owned docstrings
        # (the term author's view of the interface), not only ERGM's
        for nm in (:name, :compute, :change_stat)
            @test haskey(meta, Base.Docs.Binding(ERGMUserterms, nm))
            @test any(startswith(b, string(nm)) for (b, _) in blocks)
        end
        @test length(blocks) >= 18
        for (nm, code) in blocks
            m = Module(Symbol("DocExample_", nm))
            ok = try
                Core.eval(m, :(using ERGMUserterms))
                with_logger(NullLogger()) do
                    redirect_stdout(devnull) do
                        Core.eval(m, Meta.parseall(code; filename="docstring:$nm"))
                    end
                end
                true
            catch err
                println(stderr, "docstring example of $nm failed: ", sprint(showerror, err))
                false
            end
            @test ok
        end
    end

    @testset "README and docs pages execute without warnings" begin
        # The CI-runnable form of the site's tools/check_snippets.jl (which
        # needs a `../.snippet-env` with all 15 packages dev'ed — nothing CI
        # can build): the same extraction rules, the same execution model
        # (each file's blocks in order in one fresh module, REPL soft scope,
        # inside a temp dir, stdout swallowed), and one loophole closed —
        # every validator warns on every failure, so a `:warn` record from
        # ERGMUserterms while a page runs means an example that "ran" but was
        # false. That is a test failure here, not a line in a log.
        pkg = dirname(@__DIR__)
        pages = [joinpath(pkg, "README.md")]
        for (dir, _, files) in walkdir(joinpath(pkg, "docs", "src")), f in files
            endswith(f, ".md") && push!(pages, joinpath(dir, f))
        end
        @test length(pages) >= 7

        # Grep guards: the pre-rename module name must not come
        # back; the docs model the rng contract (no bare draws); every entry
        # page teaches the one import idiom
        all_lines = [l for p in pages for l in readlines(p)]
        @test isempty(filter(l -> occursin(r"using .*\bNetwork\b", l), all_lines))
        @test isempty(filter(l -> occursin(r"\brand\(1:", l), all_lines))
        for p in ("README.md", "docs/src/index.md", "docs/src/getting_started.md")
            @test occursin("import ERGM: name, compute, change_stat",
                           _readtext(joinpath(pkg, p)))
        end
        for p in pages
            text = _readtext(p)
            occursin(r"\b(validate_term|test_term|benchmark_term)\(", text) &&
                @test occursin("rng=Xoshiro(", text)
        end

        # --- the extractor: check_snippets.jl's rules ------------------------
        output_label = r"(?:^|\b)(?:example\s+)?output:?\**\s*$"i
        function snippets(path)
            lines = readlines(path)
            out = Tuple{Int, String, Union{Nothing, String}}[]
            i = 1
            while i <= length(lines)
                m = match(r"^\s*```(julia[a-z-]*)\s*$", lines[i])
                if m === nothing
                    i += 1
                    continue
                end
                lang = m.captures[1]
                fence = i
                body = String[]
                i += 1
                while i <= length(lines) && !occursin(r"^\s*```\s*$", lines[i])
                    push!(body, lines[i])
                    i += 1
                end
                i += 1
                code = join(body, "\n")
                skip = if lang != "julia"
                    "fenced as $lang"
                elseif fence > 1 && occursin("<!-- skip-check -->", lines[fence - 1])
                    "skip-check marker"
                elseif occursin(r"\bPkg\.(add|develop|rm|activate|instantiate|update)\(", code)
                    "Pkg install instructions"
                else
                    j = fence - 1
                    while j >= 1 && isempty(strip(lines[j]))
                        j -= 1
                    end
                    (j >= 1 && occursin(output_label, strip(lines[j]))) ? "output-only block" : nothing
                end
                push!(out, (fence, code, skip))
            end
            return out
        end
        has_parse_error(ex) = ex isa Expr &&
            (ex.head in (:error, :incomplete) || any(has_parse_error, ex.args))

        # A block that re-lists its `using` lines is what a reader pastes into
        # a fresh REPL: if it uses `Xoshiro(` it must also `using Random` (a
        # page runs sequentially under the checker, so an earlier block's
        # `using Random` would otherwise hide the omission, as it once did
        # on getting_started.md's "Complete Example")
        not_self_contained = String[]
        for p in pages, (fence, code, skip) in snippets(p)
            skip === nothing || continue
            occursin("Xoshiro(", code) && occursin(r"^\s*using\b"m, code) || continue
            occursin(r"^\s*using\b[^\n]*\bRandom\b"m, code) ||
                push!(not_self_contained, "$(relpath(p, pkg)):$fence")
        end
        @test isempty(not_self_contained)

        failures = String[]
        n_run = Ref(0)
        n_claims = Ref(0)
        # Output comments are not executed by any gate, so a page can quote
        # a number its own code does not produce (as they once did:
        # validation.md said `compute() returns: 30.0` for a 26.0 network,
        # getting_started.md `4.0` for 3.0). Every numeric
        # `compute() returns: X` / `Works on empty|complete network: X`
        # claim on a page must be a value some run ON THAT PAGE logged or
        # printed while the page executed.
        claim_re = r"(compute\(\) returns|Works on (?:empty|complete) network): (-?\d+(?:\.\d+)?)"
        for path in pages
            rel = relpath(path, pkg)
            sandbox = Module(gensym(basename(path)))
            # the checker's self-referential eval/include, so a page may
            # `include` a file (validation.md includes the package template)
            Core.eval(sandbox, :(eval(x) = Core.eval($sandbox, x)))
            Core.eval(sandbox, :(include(f) = Base.include($sandbox, f)))
            stdout_path, stdout_io = mktemp()
            logs, _ = Test.collect_test_logs(; min_level=Logging.Info) do
                cd(mktempdir()) do
                    for (fence, code, skip) in snippets(path)
                        skip === nothing || continue
                        parsed = try
                            Meta.parseall(code; filename=path)
                        catch
                            nothing
                        end
                        # unparseable blocks are output, as in the checker
                        (parsed === nothing || has_parse_error(parsed)) && continue
                        n_run[] += 1
                        stmts = parsed isa Expr && parsed.head == :toplevel ? parsed.args : Any[parsed]
                        line = fence
                        ok = true
                        for st in stmts
                            if st isa LineNumberNode
                                line = fence + st.line
                                continue
                            end
                            try
                                redirect_stdout(stdout_io) do
                                    Core.eval(sandbox, REPL.softscope(st))
                                end
                            catch err
                                push!(failures, "$rel:$line: " *
                                      first(split(sprint(showerror, err), '\n')))
                                ok = false
                                break
                            end
                        end
                        ok || break     # later blocks depend on earlier state
                    end
                end
            end
            for r in logs
                r.level >= Logging.Warn && r._module === ERGMUserterms &&
                    push!(failures, "$rel: ERGMUserterms warned while the page ran " *
                                    "(a claim on the page is false): $(r.message)")
            end
            close(stdout_io)
            produced = Set{String}()
            for text in (_readtext(stdout_path), (r.message for r in logs)...)
                for m in eachmatch(claim_re, text)
                    push!(produced, m.captures[1] * ": " * m.captures[2])
                end
            end
            for (i, line) in enumerate(readlines(path)), m in eachmatch(claim_re, line)
                n_claims[] += 1
                claim = m.captures[1] * ": " * m.captures[2]
                claim in produced ||
                    push!(failures, "$rel:$i quotes `$claim` but no run on the page " *
                                    "produced that value")
            end
        end
        @test n_run[] >= 50
        @test n_claims[] >= 6      # the four compute() and two Works-on claims
        foreach(f -> println(stderr, "docs example failed: ", f), failures)
        @test isempty(failures)
    end
    # ------------------------------------------------------------------
    # Golden fixture against statnet
    # ------------------------------------------------------------------
    @testset "Golden: bundled terms and the harness match statnet" begin
        # `load_golden` throws on a fixture without a [provenance] block or
        # whose generating script is gone — that is the gate. Every value
        # below was emitted by ergm 4.12.0 (test/fixtures/r/userterms_examples.R).
        g = load_golden(joinpath(@__DIR__, "fixtures", "userterms_examples.toml"))
        @test g.provenance["ergm_version"] == "4.12.0"
        @test isfile(g.script_path)
        v = g.values
        n = Int(v["n"])

        # Rebuild R's network from the TOML alone: arcs, the three vertex
        # attributes and the five stored weights. Every OTHER arc — the
        # reverse arcs (3,1), (4,2), (5,2) included — has no weight attribute
        # and must count at WeightedEdges' default, which is W's 1.0 in R.
        function golden_net(directed::Bool)
            net = network(n; directed=directed)
            tails, heads = directed ? (v["tails"], v["heads"]) :
                                      (v["und_tails"], v["und_heads"])
            for (i, j) in zip(tails, heads)
                add_edge!(net, i, j)
            end
            set_vertex_attribute!(net, :a, Dict(k => Float64(v["a"][k]) for k in 1:n))
            set_vertex_attribute!(net, :b, Dict(k => Float64(v["b"][k]) for k in 1:n))
            set_vertex_attribute!(net, :group, Dict(k => String(v["group"][k]) for k in 1:n))
            for (i, j) in zip(v["weighted_tails"], v["weighted_heads"])
                set_edge_attribute!(net, :weight, i, j, Float64(v["weight"]))
            end
            return net
        end
        netd, netu = golden_net(true), golden_net(false)
        @test ne(netd) == 11 && ne(netu) == 8

        # The dyadic covariate is the ASYMMETRIC cov[i,j] = 2i + j on both
        # networks: R's undirected block got it symmetrised to its (min,max)
        # entry, which is exactly the canonical key `DyadCovTerm` reads on an
        # undirected network — so agreement pins that branch.
        cov = [Float64(2i + j) for i in 1:n, j in 1:n]
        terms = [("example",         ExampleTerm()),            # nodecov("id")
                 ("template",        TemplateTerm(1.0)),        # edges
                 ("weightededges",   WeightedEdges()),          # edgecov(W)
                 ("dyadcov",         DyadCovTerm(cov)),         # edgecov(cov)
                 ("interaction",     InteractionTerm(:a, :b)),  # edgecov(M)
                 ("recip_homophily", ReciprocatedHomophily(:group))]  # mutual(same="group")
        und_terms = terms[1:5]        # mutual is directed-only, as in R

        # (a) full statistics
        @test check_golden(g, "summary_directed", [compute(t, netd) for (_, t) in terms]) ||
              println(golden_report(g, "summary_directed", [compute(t, netd) for (_, t) in terms]))
        @test check_golden(g, "summary_undirected", [compute(t, netu) for (_, t) in und_terms]) ||
              println(golden_report(g, "summary_undirected", [compute(t, netu) for (_, t) in und_terms]))

        # (b)/(c) every add-direction change statistic, in the fixture's
        # documented order (row-major over ordered dyads i != j; over i < j
        # undirected), from the term's own `change_stat` AND from the
        # harness's brute-force toggle-and-recompute reference — so the
        # reference `validate_term` compares user terms against is itself
        # R's number on every dyad, not a self-consistency check.
        dyads_d = [(i, j) for i in 1:n for j in 1:n if i != j]
        dyads_u = [(i, j) for i in 1:n for j in 1:n if i < j]
        brute = ERGMUserterms._brute_change_stat
        for (key, t) in terms
            k = "change_$(key)_directed"
            @test check_golden(g, k, [change_stat(t, netd, i, j) for (i, j) in dyads_d]) ||
                  println(golden_report(g, k, [change_stat(t, netd, i, j) for (i, j) in dyads_d]))
            @test check_golden(g, k, [brute(t, netd, i, j) for (i, j) in dyads_d]) ||
                  println(golden_report(g, k, [brute(t, netd, i, j) for (i, j) in dyads_d]))
        end
        for (key, t) in und_terms
            k = "change_$(key)_undirected"
            @test check_golden(g, k, [change_stat(t, netu, i, j) for (i, j) in dyads_u])
            # ... and asked in the (j, i) order, which is where the (min,max)
            # canonical-key branches of WeightedEdges/DyadCovTerm live
            @test check_golden(g, k, [change_stat(t, netu, j, i) for (i, j) in dyads_u]) ||
                  println(golden_report(g, k, [change_stat(t, netu, j, i) for (i, j) in dyads_u]))
            @test check_golden(g, k, [brute(t, netu, i, j) for (i, j) in dyads_u])
        end
        # The brute-force reference restored everything it toggled
        @test ne(netd) == 11 && ne(netu) == 8
        @test get_edge_attribute(netd, :weight, 1, 3) == 2.5
        @test get_edge_attribute(netd, :weight, 3, 1) === nothing
        @test get_edge_attribute(netu, :weight, 3, 1) == 2.5

        # (d) WeightedEdges as a FITTED model is edgecov(W) — by MPLE and by
        # MCMLE alike. The statistics above cannot see a term that reads the
        # edge attribute live: ERGM's samplers remove and re-add edges, and
        # rem_edge! deletes the stored weight, so before the snapshot
        # (`ERGM.Extension.materialize(::WeightedEdges, net)`) every weight was lost
        # as the chain ran. The model is dyad-independent, so R's value is
        # the exact MLE and an MCMLE fit must land on it.
        for (d, net) in (("directed", netd), ("undirected", netu))
            mple = fit_ergm(net, [WeightedEdges()]; method=:mple)
            @test check_golden(g, "estimate_weightededges_$d", coef(mple)[1]) ||
                  println(golden_report(g, "estimate_weightededges_$d", coef(mple)[1]))
            @test check_golden(g, "exact_estimate_weightededges_$d", coef(mple)[1])
            @test check_golden(g, "se_weightededges_$d", stderror(mple)[1])
            @test check_golden(g, "exact_se_weightededges_$d", stderror(mple)[1])
            @test check_golden(g, "loglik_weightededges_$d", loglikelihood(mple))
            mcmle = fit_ergm(net, [WeightedEdges()]; method=:mcmle, rng=Xoshiro(7))
            @test mcmle.converged
            @test check_golden(g, "estimate_weightededges_$d", coef(mcmle)[1]) ||
                  println(golden_report(g, "estimate_weightededges_$d", coef(mcmle)[1]))
            @test check_golden(g, "loglik_weightededges_$d", loglikelihood(mcmle))
            # The MCMLE standard error is a Monte-Carlo estimate from draws of
            # the model: on R's value to Monte-Carlo accuracy (measured 3 %)
            @test isapprox(stderror(mcmle)[1], v["exact_se_weightededges_$d"]; rtol=0.15)
            # the fits left the caller's weights alone
            @test length(get_edge_attribute(net, :weight)) == 5
        end

        # The R fixture's coefficient labels are recorded but NOT claimed:
        # user terms keep their own names (a term author chooses the label,
        # and `nodecov.id` would be a lie about what ExampleTerm reads)
        r_names = String.(v["r_names_directed"])
        @test r_names == ["nodecov.id", "edges", "edgecov.W", "edgecov.cov", "edgecov.M", "mutual.group"]
        @test [name(t) for (_, t) in terms] ==
              ["example", "template.1.0", "weightedges.weight", "dyadcov", "interact.a.b", "recip_homophily.group"]
        @test isempty(intersect(r_names, [name(t) for (_, t) in terms]))
    end

    # ------------------------------------------------------------------
    # Allocation pins
    # ------------------------------------------------------------------
    @testset "Allocation regressions" begin
        # The bundled terms are what third-party authors copy, so their
        # per-dyad cost is pinned at 0 B — `compute` too — and so is
        # inference: `change_stat`/`compute` must infer `Float64`. Development
        # versions shipped `WeightedEdges` and `InteractionTerm` inferring `Any` (values
        # pulled from the untyped `Dict{…,Any}` attribute stores with `get`),
        # which poisoned ERGM's statically typed change-statistic tuple
        # (`Tuple{Float64, Any}`) for every model containing them and boxed on
        # every MPLE row and MH step. They now read per element with
        # `get_vertex_attribute(net, attr, v)` / `get_edge_attribute(net,
        # attr, i, j)` (0 B; the whole-Dict getters allocate an empty default
        # Dict on every call) and assert `::Float64`. The typed snapshot at
        # construction (`vertex_attribute_vector(net, attr, Float64)`,
        # term_interface.md) is still the faster pattern for a hot path — a
        # vector index instead of two hash lookups — but a copied template
        # must never infer `Any` again.
        n = 60
        net = random_net(n, 300; directed=true, seed=77)
        set_vertex_attribute!(net, :a, Dict(k => Float64(k) for k in 1:n))
        set_vertex_attribute!(net, :b, Dict(k => Float64(n + 1 - k) for k in 1:n))
        set_vertex_attribute!(net, :group, Dict(k => isodd(k) ? "a" : "b" for k in 1:n))
        for e in first(collect(edges(net)), 20)
            set_edge_attribute!(net, :weight, src(e), dst(e), 2.5)
        end
        cov = [Float64(2i + j) for i in 1:n, j in 1:n]
        rng = Xoshiro(78)
        dyads = [ERGMUserterms._random_dyad(rng, n) for _ in 1:40]

        # Worst case over the dyad sample, after a warm-up call per dyad
        function worst_change_alloc(t, net, dyads)
            worst = 0
            for (i, j) in dyads
                change_stat(t, net, i, j)
                worst = max(worst, @allocated change_stat(t, net, i, j))
            end
            return worst
        end
        function compute_alloc(t, net)
            compute(t, net)
            return @allocated compute(t, net)
        end

        und = random_net(n, 150; directed=false, seed=76)
        set_vertex_attribute!(und, :a, Dict(k => Float64(k) for k in 1:n))
        set_vertex_attribute!(und, :b, Dict(k => Float64(n + 1 - k) for k in 1:n))
        for e in first(collect(edges(und)), 20)
            set_edge_attribute!(und, :weight, src(e), dst(e), 2.5)
        end
        bundled = (ExampleTerm(), TemplateTerm(1.5), WeightedEdges(), DyadCovTerm(cov),
                   InteractionTerm(:a, :b), ReciprocatedHomophily(:group))
        for t in bundled
            @test worst_change_alloc(t, net, dyads) == 0
            @test compute_alloc(t, net) == 0
            # Inference closes on Float64 on both network kinds (the
            # undirected one exercises the canonical-key branches)
            for g in (net, und)
                @test Base.return_types(change_stat, (typeof(t), typeof(g), Int, Int)) == [Float64]
                @test Base.return_types(compute, (typeof(t), typeof(g))) == [Float64]
            end
        end
        for t in bundled[1:5]      # ReciprocatedHomophily is directed-only
            @test worst_change_alloc(t, und, dyads) == 0
            @test compute_alloc(t, und) == 0
        end
        # The form a model holds `WeightedEdges` in (the weight-matrix
        # snapshot ERGM's samplers and estimators evaluate) is pinned the same
        for g in (net, und)
            held = ERGM.Extension.materialize(WeightedEdges(), g)
            @test worst_change_alloc(held, g, dyads) == 0
            @test compute_alloc(held, g) == 0
            @test Base.return_types(change_stat, (typeof(held), typeof(g), Int, Int)) == [Float64]
            @test Base.return_types(compute, (typeof(held), typeof(g))) == [Float64]
        end

        # The harness's own building blocks stay O(n²) at worst:
        # `_term_fingerprint` is one `compute` plus exactly one `change_stat`
        # per dyad, and `_other_dyads` caps the dyads it returns at 250
        # however large the network (30 vertices have 870 ordered dyads).
        small = random_net(8, 12; directed=true, seed=79)
        t = RecordingTerm()
        ERGMUserterms._term_fingerprint(t, small)
        @test length(t.visited) == 8 * 7
        @test allunique(t.visited)
        @test length(ERGMUserterms._other_dyads(Xoshiro(1), small, 1, 2)) == 8 * 7 - 1
        big = network(30; directed=true)
        for k in 1:29
            add_edge!(big, k, k + 1)
        end
        @test length(ERGMUserterms._other_dyads(Xoshiro(1), big, 1, 2)) <= 250
        @test length(ERGMUserterms._other_dyads(Xoshiro(1), big, 1, 2)) == 250
        @test (1, 2) ∉ ERGMUserterms._other_dyads(Xoshiro(1), big, 1, 2)
    end

    # ------------------------------------------------------------------
    # The harness's power, g(∅), trusted declarations and edge attributes
    # under MCMC
    # ------------------------------------------------------------------
    @testset "validate_term checks every dyad by default (power)" begin
        # A mutual-like term whose change statistic is wrong whenever i > j —
        # on half of all reciprocated dyads. With the old default (10 uniform
        # random dyads) it PASSED on 136 of 200 random 20-vertex networks at
        # density 0.1.
        function sparse_net(seed)
            rng = Xoshiro(seed)
            net = network(20; directed=true)
            for i in 1:20, j in 1:20
                i != j && rand(rng) < 0.1 && add_edge!(net, i, j)
            end
            return net
        end
        nets = [sparse_net(s) for s in 1:200]
        # (the bug is observable on every one of them: any arc i -> j with
        # i > j has a reverse dyad whose change statistic is wrong)
        passes = [validate_term(HalfWrongMutual(), net; verbose=false, traits=false)
                  for net in nets]
        @test !any(passes)
        # The opt-in sample is stratified by dyad state, so even 10 dyads
        # reach the reciprocated ones (measured: 1 pass in 200; 136 before)
        sampled = [validate_term(HalfWrongMutual(), nets[s]; verbose=false, traits=false,
                                 n_tests=10, rng=Xoshiro(1000 + s)) for s in 1:200]
        @test count(sampled) <= 10
        # The verdict line says how much was checked
        logs = Test.collect_test_logs(() -> validate_term(ExampleTerm(), nets[1]; traits=false))[1]
        @test any(occursin("consistent with compute() on all 380 dyads", l.message) for l in logs)
        logs = Test.collect_test_logs(() -> validate_term(ExampleTerm(), nets[1]; traits=false,
                                                          n_tests=25, rng=Xoshiro(1)))[1]
        @test any(occursin("on 25 of 380 dyads (a stratified sample", l.message) for l in logs)

        # Stratification: every non-empty stratum is represented, small
        # strata are exhausted first, no dyad twice
        net = nets[1]
        all_dyads = ERGMUserterms._all_dyads(net)
        strata = Dict{Any, Int}()
        for d in all_dyads
            k = ERGMUserterms._dyad_stratum(net, d...)
            strata[k] = get(strata, k, 0) + 1
        end
        pick = ERGMUserterms._stratified_dyads(Xoshiro(3), net, all_dyads, 40)
        @test length(pick) == 40 && allunique(pick)
        got = Dict{Any, Int}()
        for d in pick
            k = ERGMUserterms._dyad_stratum(net, d...)
            got[k] = get(got, k, 0) + 1
        end
        @test keys(got) == keys(strata)
        share = 40 ÷ length(strata)
        @test all(got[k] >= min(strata[k], share) for k in keys(strata))
        @test ERGMUserterms._stratified_dyads(Xoshiro(3), net, all_dyads, 10_000) == all_dyads
        # A network above the exhaustive limit gets a stratified 5000, and says so
        big = random_net(80, 300; seed=4)
        dy, words = ERGMUserterms._dyads_to_check(Xoshiro(1), big, nothing)
        @test length(dy) == 5000 && allunique(dy)
        @test occursin("5000 of 6320 dyads (a stratified sample; the network has more than 5000 dyads)", words)
    end

    @testset "g(∅) and the path from the empty network are checked" begin
        flo = load_dataset(:florentine_marriage)
        # The old guide's advice — `ne(net) == 0 && return 0.0` — gives an
        # isolates count that is wrong on the empty network (g(∅) = n, not 0).
        # No per-dyad check on the observed network sees it; it used to pass.
        @test !validate_term(GuideIsolates(), flo; verbose=false)
        @test_logs (:warn, r"compute\(\) on the empty network \(g\(∅\) = 0\.0\) plus the change statistics summed while the 20 edges are added one at a time gives -15\.0, but compute\(\) on the network gives 1\.0") match_mode=:any validate_term(GuideIsolates(), flo; traits=false)
        # The same statistic without the special case passes, with g(∅) = 16
        @test validate_term(PlainIsolates(), flo; verbose=false)
        @test_logs (:info, r"✓ g\(∅\) = 16\.0 and the summed change statistics from the empty network reproduce compute\(\)") match_mode=:any validate_term(PlainIsolates(), flo; traits=false)
        @test compute(PlainIsolates(), network(16; directed=false)) == 16.0
        @test compute(PlainIsolates(), flo) == 1.0
        # ... and a term with g(∅) ≠ 0 is fitted with the right likelihood:
        # ERGM's exact log-normaliser includes θ'g(∅). The number of non-ties
        # is a reparametrisation of the edge count (nonties = 120 − edges), so
        # its MCMLE log-likelihood — exact reference plus bridge — must be the
        # Bernoulli log-likelihood at the observed density, computed here by
        # brute force over the 120 dyads; without θ'g(∅) it was off by
        # θ̂ · 120 (about +194)
        own = fit_ergm(flo, [NonTies()]; method=:mcmle, rng=Xoshiro(4))
        p̂ = ne(flo) / 120
        @test coef(own)[1] ≈ -log(p̂ / (1 - p̂)) atol = 1e-8
        bernoulli = ne(flo) * log(p̂) + (120 - ne(flo)) * log(1 - p̂)
        @test loglikelihood(own) ≈ bernoulli atol = 1e-8
        @test loglikelihood(fit_ergm(flo, [NonTies()]; method=:mple)) ≈ bernoulli atol = 1e-8
        # A wrong change statistic away from the observed network also
        # fails the path check (here: wrong only on the empty network)
        @test !validate_term(WrongWhenEmpty(), flo; verbose=false, traits=false)
        @test change_stat_check(WrongWhenEmpty(), flo; verbose=false)   # per-dyad check cannot see it
        # Edge attributes are carried along the path (exogenous dyadic data)
        w = deepcopy(flo)
        for (k, e) in enumerate(collect(edges(flo)))
            set_edge_attribute!(w, :weight, src(e), dst(e), 0.5 + k)
        end
        @test ERGMUserterms._path_consistent(WeightedEdges(), w, 1e-10, false)
        @test length(get_edge_attribute(w, :weight)) == 20
    end

    @testset "A mis-declared dyad-independence always fails" begin
        # `is_exact`, the exact log-normaliser and "MPLE = MLE" trust the
        # declaration. The check used to toggle around 4 random source dyads
        # and caught a mutual-type term on 82 of 100 seeds.
        sam = load_dataset(:sampson)
        @test all(!validate_traits(LyingDependenceTerm(), sam; verbose=false, rng=Xoshiro(s))
                  for s in 1:100)
        @test_logs (:warn, r"declares is_dyad_dependent = false, but its change statistic at dyad \(\d+,\d+\) is .* on the network but .* on the empty network: it IS dyad-dependent") match_mode=:any validate_traits(LyingDependenceTerm(), sam; rng=Xoshiro(1))
        @test !validate_term(LyingDependenceTerm(), sam; verbose=false)
        # A dependence invisible on the empty and the complete network —
        # the change statistic reads ONE distant dyad — is found by the
        # exhaustive pairwise toggling
        @test !validate_traits(ReadsOneDyad(), sam; verbose=false, rng=Xoshiro(1))
        @test_logs (:warn, r"moved from .* when dyad \(17,18\) was toggled") match_mode=:any validate_traits(ReadsOneDyad(), sam; rng=Xoshiro(1))
        # On a network too big for all pairs, the reverse arc and the
        # incident dyads are toggled first: a mutual-type lie is still found
        # at every seed
        big = random_net(40, 150; seed=8)
        @test all(!validate_traits(LyingDependenceTerm(), big; verbose=false, rng=Xoshiro(s))
                  for s in 1:20)
        near = ERGMUserterms._other_dyads(Xoshiro(1), big, 3, 9; cap=250)
        @test (9, 3) in near
        @test count(d -> any(in((3, 9)), d), near) == 4 * 38 + 1
        # Truthful declarations keep passing, as does the declared-dependent
        # default (which asserts nothing)
        @test validate_traits(ExampleTerm(), sam; verbose=false, rng=Xoshiro(1))
        @test validate_traits(SharedNeighborTerm(), sam; verbose=false, rng=Xoshiro(1))
    end

    @testset "WeightedEdges is edgecov(W) under MCMC" begin
        # A 20-vertex digraph with 49 weighted arcs (weights in [0.5, 3])
        rng = Xoshiro(42)
        n = 20
        net = network(n; directed=true)
        W = ones(n, n)
        for i in 1:n, j in 1:n
            i == j && continue
            if rand(rng) < 0.15
                add_edge!(net, i, j)
                w = 0.5 + 2.5rand(rng)
                set_edge_attribute!(net, :weight, i, j, w)
                W[i, j] = w
            end
        end
        @test ne(net) == 49
        t = WeightedEdges()
        model = ERGMModel(ERGMFormula([t]), net)
        held = model.formula.terms[1]
        @test held isa ERGMUserterms.MaterializedWeightedEdges
        @test held.weights == W
        @test name(held) == name(t) == "weightedges.weight"
        @test !ERGM.is_dyad_dependent(held)
        @test ERGM.Extension.materialize(held, net) === held
        @test compute(held, net) == compute(t, net) == compute(EdgeCov(W), net)

        # The sampler deletes every stored weight (that is what rem_edge!
        # does) — and the model does not care: draw for draw, the chain is
        # the chain of EdgeCov(W)
        θ = [-1.0]
        sims = sample_networks(model, θ; n_sim=500, rng=Xoshiro(1), n_chains=1)
        ref = sample_networks(ERGMModel(ERGMFormula([EdgeCov(W)]), net), θ;
                              n_sim=500, rng=Xoshiro(1), n_chains=1)
        edge_set(s) = Set((Int(src(e)), Int(dst(e))) for e in edges(s))
        @test all(edge_set(a) == edge_set(b) for (a, b) in zip(sims, ref))
        @test sum(length(get_edge_attribute(s, :weight)) for s in sims) < 49 * 500 ÷ 2
        @test [compute(held, s) for s in sims] == [compute(EdgeCov(W), s) for s in sims]

        # Power: with weight 3 on the observed arcs the edgecov(W) model and
        # the all-default model the old live-reading term decayed to are far
        # apart, and the chain sits on the former
        heavy = deepcopy(net)
        W3 = ones(n, n)
        for e in collect(edges(net))
            set_edge_attribute!(heavy, :weight, src(e), dst(e), 3.0)
            W3[src(e), dst(e)] = 3.0
        end
        expect(M) = sum(M[i, j] / (1 + exp(M[i, j])) for i in 1:n, j in 1:n if i != j)
        variance(M) = sum(M[i, j]^2 * exp(M[i, j]) / (1 + exp(M[i, j]))^2
                          for i in 1:n, j in 1:n if i != j)
        target, decayed = expect(W3), expect(ones(n, n))
        @test decayed - target > 6                       # 95.99 vs 102.20
        hmodel = ERGMModel(ERGMFormula([t]), heavy)
        hsims = sample_networks(hmodel, θ; n_sim=2000, rng=Xoshiro(2), n_chains=1)
        vals = [compute(hmodel.formula.terms[1], s) for s in hsims]
        # lag-1 autocorrelation of the draws is ≈ 0.6: ESS ≈ n/4
        mcse = sqrt(variance(W3) / (2000 / 4))
        @test abs(sum(vals) / 2000 - target) < 4 * mcse
        @test abs(sum(vals) / 2000 - decayed) > 8 * mcse

        # Estimation: dyad-independent, so the exact MLE is the logistic
        # MLE (= R's edgecov(W) fit). Brute force: the one-parameter score
        # equation Σ W_ij y_ij = Σ W_ij expit(θ W_ij) over the 380 dyads,
        # solved by bisection (the score is decreasing in θ). The MCMLE used
        # to return −1.266 with converged = true.
        score(θ) = sum(W[i, j] * ((has_edge(net, i, j) ? 1.0 : 0.0) -
                                  1 / (1 + exp(-θ * W[i, j])))
                       for i in 1:n, j in 1:n if i != j)
        lo, hi = -5.0, 5.0
        for _ in 1:200
            mid = (lo + hi) / 2
            score(mid) > 0 ? (lo = mid) : (hi = mid)
        end
        θ̂ = (lo + hi) / 2
        @test abs(score(θ̂)) < 1e-9
        mple = fit_ergm(net, [t]; method=:mple)
        @test coef(mple)[1] ≈ θ̂ atol=1e-8
        @test coef(mple) ≈ coef(fit_ergm(net, [EdgeCov(W)]; method=:mple)) atol=1e-12
        mcmle = fit_ergm(net, [t]; method=:mcmle, rng=Xoshiro(3))
        @test mcmle.converged
        @test coef(mcmle)[1] ≈ θ̂ atol=1e-6
        @test is_exact(mple)
        # simulate_ergm and gof go through the same model
        draws = simulate_ergm(mple; n_sim=3, rng=Xoshiro(5))
        @test length(draws) == 3
        @test length(get_edge_attribute(net, :weight)) == 49       # caller's network untouched

        # An undirected network: the snapshot is symmetric
        und = network(5; directed=false)
        add_edge!(und, 1, 2); add_edge!(und, 2, 3)
        set_edge_attribute!(und, :weight, 2, 1, 4.0)
        mu = ERGM.Extension.materialize(t, und)
        @test mu.weights[1, 2] == mu.weights[2, 1] == 4.0
        @test change_stat(mu, und, 2, 1) == change_stat(t, und, 2, 1) == 4.0
        @test compute(mu, und) == compute(t, und) == 5.0
        # Non-numeric weights and a snapshot used on another network are
        # ArgumentErrors
        bad = deepcopy(und); set_edge_attribute!(bad, :weight, 2, 3, "heavy")
        @test_throws ArgumentError ERGMModel(ERGMFormula([t]), bad)
        @test_throws ArgumentError compute(mu, network(6; directed=false))
        @test_throws ArgumentError change_stat(mu, network(6; directed=false), 1, 2)

        # The harness now sees the trap: a term that reads the attribute
        # live passes every numeric check (the harness restores attributes
        # around its toggles) and FAILS the sampler check; the bundled term
        # passes
        @test validate_term(t, net; verbose=false)
        @test change_stat_check(LiveWeight(), net; verbose=false)
        @test !validate_term(LiveWeight(), net; verbose=false)
        @test_logs (:warn, r"statistic changes when the network's edges are toggled off and on again without their edge attributes being restored.*ERGM\.Extension\.materialize") match_mode=:any validate_traits(LiveWeight(), net; rng=Xoshiro(1))
        # ... a ✓ only when the check ran: no edge attributes, no line
        plain = random_net(8, 12; seed=3)
        logs = Test.collect_test_logs(() -> validate_traits(ExampleTerm(), plain; rng=Xoshiro(1)))[1]
        @test !any(occursin("sampler's toggles", l.message) for l in logs)
        logs = Test.collect_test_logs(() -> validate_traits(t, net; rng=Xoshiro(1)))[1]
        @test any(occursin("✓ statistic survives the sampler's toggles", l.message) for l in logs)
    end

    @testset "DyadCovTerm refuses a wrongly sized matrix" begin
        net = random_net(16, 40; seed=6)
        small = DyadCovTerm(ones(4, 4))
        # was: entries outside the matrix counted 0, and the model fitted to -Inf
        @test_throws ArgumentError compute(small, net)
        @test_throws ArgumentError change_stat(small, net, 1, 2)
        err = try ERGMModel(ERGMFormula([Edges(), small]), net) catch e; e end
        @test err isa ArgumentError
        @test occursin("term 'dyadcov' has a 4×4 covariate matrix but the network has 16 vertices", err.msg)
        @test_throws ArgumentError fit_ergm(net, [Edges(), small])
        @test !validate_term(small, net; verbose=false)
        @test_throws ArgumentError DyadCovTerm(ones(3, 4))
        # the right size is unchanged, and an integer matrix is accepted
        ok = DyadCovTerm([2i + j for i in 1:16, j in 1:16])
        @test ok.covariate isa Matrix{Float64}
        @test validate_term(ok, net; verbose=false)
        @test ERGMModel(ERGMFormula([Edges(), ok]), net).formula.terms[2] === ok
    end

    @testset "README limitations, CHANGELOG and the executed README term agree" begin
        root = dirname(@__DIR__)
        readme = _readtext(joinpath(root, "README.md"))
        index = _readtext(joinpath(root, "docs", "src", "index.md"))
        changelog = _readtext(joinpath(root, "CHANGELOG.md"))
        section(text, heading) = (m = match(Regex("\\n#{2,3} " * heading * "[^\\n]*\\n(.*?)(?=\\n## )", "s"), text);
                                  m === nothing ? "" : m.captures[1])
        not_impl = section(readme, "Not implemented")
        not_impl_docs = section(index, "Not implemented")
        known = section(changelog, "Known limitations")
        @test !isempty(not_impl) && !isempty(not_impl_docs) && !isempty(known)
        # the same items in all three (normalising line breaks)
        squash(t) = replace(t, r"\s+" => " ")
        for item in ("mindegree", "5000 dyads", "about 500 dyads",
                     "refuse", "ArgumentError", "is_dyad_dependent",
                     "identifiability, non-degeneracy or model fit",
                     "configuration")
            @test occursin(item, squash(not_impl))
            @test occursin(item, squash(not_impl_docs))
            @test occursin(item, squash(known))
        end
        # The README's own example term — an UNANNOTATED change_stat(…, i, j),
        # which ERGM 0.1's `::Int`-typed fallbacks made ambiguous inside
        # fit_ergm — is an executed block, not a skip-check sketch (the
        # docs-execution testset above runs it through fit_ergm)
        m = match(r"([^\n]*)\n```julia\n(?:(?!```).)*?change_stat\(::IdSum, net, i, j\)(?:(?!```).)*?fit_ergm\(flo, \[Edges\(\), IdSum\(\)\]\)"s, readme)
        @test m !== nothing
        @test m !== nothing && !occursin("skip-check", m.captures[1])
        @test !occursin("docs-stable", readme) && !occursin("/stable/", readme)
        @test !occursin("root workspace project", readme)
        @test occursin("statistical-network-analysis-with-julia.github.io/citing/", readme)
    end

    @testset "No private cross-package reach-ins" begin
        # ERGM.jl's building blocks come from its extension API,
        # `ERGM.Extension` (semver-covered): the source holds no `ERGM._x` /
        # `NetworkCore._x` reach-in, and the snapshot methods of
        # `WeightedEdges` and `DyadCovTerm` extend the extension API's
        # `materialize`, the function `ERGMModel` calls
        src = _readtext(joinpath(@__DIR__, "..", "src", "ERGMUserterms.jl"))
        @test isempty(collect(eachmatch(r"\b(?:ERGM|NetworkCore)\._\w+", src)))
        ext = Set(Symbol(m.captures[1]) for m in eachmatch(r"\bERGM\.Extension\.(\w+)", src))
        for n in (:materialize, :validate_formula, :n_observed_dyads, :require_supported_network)
            @test n in ext
        end
        @test all(n -> Base.isexported(ERGM.Extension, n), ext)
        @test hasmethod(ERGM.Extension.materialize, Tuple{WeightedEdges, Any})
        @test hasmethod(ERGM.Extension.materialize, Tuple{DyadCovTerm, Any})
        @test parentmodule(ERGM.Extension.materialize) === ERGM.Extension
    end

    @testset "Aqua" begin
        Aqua.test_all(ERGMUserterms; ambiguities=false)
        @test isempty(Test.detect_ambiguities(ERGMUserterms))
    end

    # ------------------------------------------------------------------
    # Release engineering: the CI clone lists cannot
    # drift from [sources] because they are DERIVED from it — and this pins
    # that the derivation the workflows run yields exactly the [sources] keys
    # ------------------------------------------------------------------
    @testset "Workflows reconstruct the ecosystem layout from [sources]" begin
        pkgdir = dirname(@__DIR__)
        pkg = pkgdir
        PKG = "ERGMUserterms"
        EXPECTED = Set(["ERGMUserterms.jl", "ERGM.jl", "NetworkCore.jl"])
        layout = layout_siblings(pkgdir, EXPECTED, PKG)
        layout.in_layout || @info "Workflow layout step not run: none of the sibling " *
            "checkouts $(join(layout.siblings, ", ")) is beside $(dirname(pkgdir)) (a " *
            "lone checkout or a registry install, where [sources] is not used)."
        # The predicate itself: no sibling → skip; any sibling → run (and a
        # missing one then fails the exact-set assertion)
        mktempdir() do root
            fake = joinpath(root, "$PKG.jl")
            mkpath(fake)
            @test !layout_siblings(fake, EXPECTED, PKG).in_layout
            mkpath(joinpath(root, "ERGM.jl"))
            @test !layout_siblings(fake, EXPECTED, PKG).in_layout      # a bare directory is no checkout
            touch(joinpath(root, "ERGM.jl", "Project.toml"))
            l = layout_siblings(fake, EXPECTED, PKG)
            @test l.in_layout && l.present == ["ERGM.jl"] && l.siblings == ["ERGM.jl", "NetworkCore.jl"]
        end
        sources = Set(keys(TOML.parsefile(joinpath(pkg, "Project.toml"))["sources"]))
        docs_sources = setdiff(
            Set(keys(TOML.parsefile(joinpath(pkg, "docs", "Project.toml"))["sources"])),
            ["ERGMUserterms"])
        @test sources == Set(["NetworkCore", "ERGM"])
        @test docs_sources == sources
        for wf in ("CI.yml", "Documentation.yml")
            yml = _readtext(joinpath(pkgdir, ".github", "workflows", wf))
            @test !occursin(r"for pkg in", yml)              # no hand-kept clone list
            @test !occursin("checkout_sources.jl", yml)
            @test occursin("path: $PKG.jl\n", yml)
            step = match(r"\n      - name: Reconstruct the ecosystem layout from \[sources\]\n        shell: julia[^\n]*\n        run: \|\n((?:          [^\n]*\n)+)", yml)
            @test step !== nothing
            step === nothing && continue
            @test first(findfirst("setup-julia", yml)) < step.offset
            # Run the workflow's own step without cloning: in the layout this
            # suite runs in, it must find exactly the siblings [sources] names.
            if !layout.in_layout
                @test_skip layout.in_layout
                continue
            end
            script = replace(step.captures[1], r"^          "m => "")
            out = mktemp() do path, io
                write(io, script); close(io)
                withenv("GITHUB_WORKSPACE" => dirname(pkgdir),
                        "GITHUB_REPOSITORY" => "statistical-network-analysis-with-Julia/$PKG.jl",
                        "LAYOUT_CHECK_ONLY" => "true", "GITHUB_STEP_SUMMARY" => nothing,
                        # `Pkg.test` runs this suite with a sandbox load path
                        # that hides the stdlibs the step loads (`TOML`)
                        "JULIA_LOAD_PATH" => nothing, "JULIA_PROJECT" => nothing) do
                    read(`$(Base.julia_cmd()) --startup-file=no $path`, String)
                end
            end
            @test Set(m.captures[1] for m in eachmatch(r"^\| (\S+\.jl) \|"m, out)) == EXPECTED
        end

        # The regression gates and the template's tests run in CI
        ci = _readtext(joinpath(pkg, ".github", "workflows", "CI.yml"))
        @test occursin("julia --project=benchmark benchmark/regression_tests.jl", ci)
        @test occursin("julia --project=examples/MyTermPackage", ci)
        @test isfile(joinpath(pkg, "benchmark", "regression_tests.jl"))
        @test isfile(joinpath(pkg, "benchmark", "benchmarks.jl"))
        bench_sources = TOML.parsefile(joinpath(pkg, "benchmark", "Project.toml"))["sources"]
        @test bench_sources["ERGMUserterms"]["path"] == ".."
        # (the sibling checkouts exist only in the layout)
        layout.in_layout &&
            @test all(isfile(joinpath(pkg, "benchmark", bench_sources[k]["path"], "Project.toml"))
                      for k in keys(bench_sources))
    end
end
