using Test
using JuMP
using GLPK
using TriObjectiveMatheuristic

@testset "Pareto filtering" begin
    X = [[1], [2], [3], [4]]
    Y = [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0],
         [1.0, 2.0, 1.0], [1.0, 1.0, 1.0]]
    _, filtered = pareto_filter(X, Y)
    @test length(filtered) == 2
    @test [2.0, 1.0, 1.0] in filtered
    @test [1.0, 2.0, 1.0] in filtered
end

@testset "Fractionality dispatch" begin
    integral = [[0.0, 1.0], [1.0, 0.0]]
    fractional = [[0.5, 0.0], [1.0, 0.5]]
    at_threshold = [[0.5, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]]
    @test first(select_algorithm(integral)) == :FPplus
    @test first(select_algorithm(fractional)) == :FPGPR
    @test first(select_algorithm(at_threshold)) == :FPplus
end

@testset "Context, feasibility and LP projection" begin
    B = [1.0 1.0]
    rhs = [1.0]
    senses = ["l"]
    C = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    ctx = build_context(B, rhs, senses, C)
    @test feasible_integer(ctx, [1, 0])
    @test !any(is_fixed, ctx.binary_variables)
    @test fbsearch(ctx, [0, 0]) !== nothing
    @test length(getobjval([1, 0], ctx)) == 3
    result = FPGPR(ctx, [[0.5, 0.0], [1.0, 0.5]], 0.0)
    @test result.method == :FPGPR
    @test result.new_solutions == 0
end

@testset "Degenerate path-relinking input" begin
    B = zeros(0, 1)
    ctx = build_context(B, Float64[], String[], ones(3, 1))
    X, Y, added = GPR(ctx, [[1.0]], [[1.0, 1.0, 1.0]], 0.1)
    @test length(X) == length(Y) == 1
    @test added == 0
end
