using ERGMUserterms
using ERGM
using Test

@testset "ERGMUserterms.jl" begin
    @testset "Module loading" begin
        @test @isdefined(ERGMUserterms)
    end

    @testset "AbstractUserTerm hierarchy" begin
        @test AbstractUserTerm <: ERGM.AbstractERGMTerm
    end

    @testset "ExampleTerm" begin
        t = ExampleTerm()
        @test t isa ExampleTerm
        @test t isa AbstractUserTerm
    end

    @testset "TemplateTerm" begin
        t = TemplateTerm(2.0)
        @test t isa TemplateTerm{Float64}
        @test t.param == 2.0
        @test t.attr == :none

        t2 = TemplateTerm(1; attr=:weight)
        @test t2.param == 1
        @test t2.attr == :weight
    end

    @testset "WeightedEdges" begin
        t = WeightedEdges()
        @test t isa WeightedEdges
        @test t.attr == :weight

        t2 = WeightedEdges(:strength)
        @test t2.attr == :strength
    end

    @testset "DyadCovTerm" begin
        cov = ones(5, 5)
        t = DyadCovTerm(cov)
        @test t isa DyadCovTerm
        @test size(t.covariate) == (5, 5)
    end

    @testset "InteractionTerm" begin
        t = InteractionTerm(:age, :income)
        @test t isa InteractionTerm
        @test t.attr1 == :age
        @test t.attr2 == :income
    end

    @testset "Validation utilities" begin
        @test isdefined(ERGMUserterms, :validate_term)
        @test isdefined(ERGMUserterms, :change_stat_check)
        @test isdefined(ERGMUserterms, :consistency_check)
        @test isdefined(ERGMUserterms, :test_term)
        @test isdefined(ERGMUserterms, :benchmark_term)
    end

    @testset "Documentation helpers" begin
        t = ExampleTerm()
        sig = term_signature(t)
        @test sig isa String
        @test occursin("ExampleTerm", sig)

        doc = term_documentation(t)
        @test doc isa String
        @test occursin("ExampleTerm", doc)
    end
end
