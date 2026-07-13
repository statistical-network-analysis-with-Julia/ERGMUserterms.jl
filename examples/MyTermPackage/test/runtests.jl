using ERGM
using ERGMUserterms
using Graphs
using MyTermPackage
using Networks
using Test

# A small directed network carrying the attribute the term declares
function example_net()
    net = network(8; directed=true)
    set_vertex_attribute!(net, :group,
                          Dict(v => isodd(v) ? "a" : "b" for v in 1:8))
    for (i, j) in ((1, 3), (3, 1), (2, 4), (4, 2), (1, 2), (5, 7), (6, 8))
        add_edge!(net, i, j)
    end
    return net
end

@testset "MyTermPackage" begin
    net = example_net()
    term = ReciprocatedHomophily(:group)

    @testset "Statistic" begin
        # Mutual dyads: (1,3) both "a", (2,4) both "b" -> 2
        @test compute(term, net) == 2.0
        @test name(term) == "recip_homophily.group"

        # Add-direction change statistic: completing a matching mutual dyad
        @test change_stat(term, net, 7, 5) == 1.0   # 5->7 exists, both "a"
        @test change_stat(term, net, 8, 6) == 1.0   # 6->8 exists, both "b"
        @test change_stat(term, net, 2, 1) == 0.0   # 1->2 exists, groups differ
        @test change_stat(term, net, 4, 5) == 0.0   # no reverse arc
    end

    @testset "Term is valid (interface + trait declarations)" begin
        # validate_term checks the add-direction convention AND every trait
        # declaration: the declared attribute is the one the term reads, the
        # direction requirement is enforced by ERGMModel, and the
        # supports_missing claim holds under face-value flips of masked dyads.
        @test validate_term(term, net; verbose=false)
    end

    @testset "Traits" begin
        @test ERGM.required_vertex_attributes(term) == (:group,)
        @test ERGM.required_edge_attributes(term) == ()
        @test ERGM.requires_directed(term)
        @test ERGM.is_dyad_dependent(term)
        @test Networks.supports_missing(term)
    end

    @testset "ERGM model construction" begin
        # Accepted, and fittable, alongside ERGM's built-in terms
        model = ERGMModel(ERGMFormula([Edges(), term]), net)
        @test model.formula.terms.names == ["edges", "recip_homophily.group"]

        fit = fit_ergm(net, [Edges(), term])
        @test length(coef(fit)) == 2

        # Declared vertex attribute absent -> the standard ArgumentError
        bare = network(8; directed=true)
        @test_throws ArgumentError ERGMModel(ERGMFormula([term]), bare)

        # Declared direction requirement -> rejected on an undirected network
        und = network(8; directed=false)
        set_vertex_attribute!(und, :group,
                              Dict(v => isodd(v) ? "a" : "b" for v in 1:8))
        @test_throws ArgumentError ERGMModel(ERGMFormula([term]), und)
    end

    @testset "Missing dyads are not counted at face value" begin
        masked = example_net()
        set_missing_dyad!(masked, 1, 3)     # the 1<->3 mutual pair is unobserved
        @test compute(term, masked) == 1.0  # only the 2<->4 pair remains
    end
end
