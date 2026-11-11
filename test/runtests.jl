using Test
using JuMP
using GLPK
using TriObjectiveMatheuristic

@testset "Pareto filtering" begin
    X = [[1], [2], [3], [4]]
    Y = [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0],
         [1.0, 2.0, 1.0], [1.0, 1.0, 1.0]]
    filtered_X, filtered_Y = pareto_filter(X, Y)
    @test length(filtered_X) == length(filtered_Y) == 2
    @test [2.0, 1.0, 1.0] in filtered_Y
    @test [1.0, 2.0, 1.0] in filtered_Y
    @test filtered_X == [[2], [3]]

    # Equal objective vectors retain one representative; empty input is valid.
    unique_X, unique_Y = pareto_filter([[1], [2]], [[1.0, 2.0, 3.0],
                                                     [1.0, 2.0, 3.0]])
    @test length(unique_X) == length(unique_Y) == 1
    @test pareto_filter(Vector{Vector{Int}}(), Vector{Vector{Float64}}()) ==
          (Vector{Vector{Int}}(), Vector{Vector{Float64}}())
    @test_throws DimensionMismatch pareto_filter([[1]], Vector{Vector{Float64}}())
end

@testset "Fractionality dispatch" begin
    integral = [[0.0, 1.0], [1.0, 0.0]]
    fractional = [[0.5, 0.0], [1.0, 0.5]]
    at_threshold = [[0.5, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]]
    below_threshold = [[0.5, 0.0], [1.0, 0.0], [0.0, 1.0],
                       [1.0, 1.0], [0.0, 0.0], [1.0, 0.0]]
    above_threshold = [[0.5, 0.0], [0.5, 1.0], [0.0, 1.0],
                       [1.0, 1.0], [0.0, 0.0]]

    @test select_algorithm(integral) == (:FPplus, 0.0)
    @test first(select_algorithm(fractional)) == :FPGPR
    @test first(select_algorithm(at_threshold)) == :FPplus
    @test first(select_algorithm(below_threshold)) == :FPplus
    @test first(select_algorithm(above_threshold)) == :FPGPR
    @test_throws ArgumentError select_algorithm(Vector{Vector{Float64}}())
    @test_throws DimensionMismatch select_algorithm([[0.0], [0.0, 1.0]])
end

@testset "Context, feasibility and LP projection" begin
    B = [1.0 1.0]
    rhs = [1.0]
    senses = ["l"]
    C = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    ctx = build_context(B, rhs, senses, C; rng=Random.MersenneTwister(17))
    @test feasible_integer(ctx, [1, 0])
    @test !feasible_integer(ctx, [0, 0])
    @test !any(is_fixed, ctx.binary_variables)
    @test fbsearch(ctx, [0, 0]) !== nothing
    @test length(getobjval([1, 0], ctx)) == 3
    @test_throws DimensionMismatch feasible_integer(ctx, [1])
    @test_throws ArgumentError feasible_integer(ctx, [0.5, 0.5])

    result = FPGPR(ctx, [[0.5, 0.0], [1.0, 0.5]], 0.0)
    @test result.method == :FPGPR
    @test result.new_solutions == 0
end

@testset "GPR end-to-end invariants" begin
    C = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    ctx = build_context(zeros(0, 2), Float64[], String[], C;
                        rng=Random.MersenneTwister(23))
    initial = [[0.0, 0.0], [1.0, 1.0]]
    objectives = [getobjval(x, ctx) for x in initial]
    X, Y, added = GPR(ctx, initial, objectives, 1.0)

    @test length(X) == length(Y)
    @test added >= 0
    @test all(feasible_integer(ctx, round.(Int, x)) for x in X)
    @test all(isapprox(Y[i], getobjval(X[i], ctx)) for i in eachindex(X))
    @test all(!dominated(Y[i], [Y[j] for j in eachindex(Y) if j != i])
              for i in eachindex(Y))
    @test_throws DimensionMismatch GPR(ctx, initial, objectives[1:1], 0.0)
end

@testset "FP+ and FPGPR dispatch end-to-end" begin
    C = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    ctx = build_context(zeros(0, 2), Float64[], String[], C;
                        rng=Random.MersenneTwister(31))

    # Integral LB vectors select FP+; all returned points must be feasible.
    fp_result = FPGPR(ctx, [[0.0, 0.0], [1.0, 1.0]], 0.25)
    @test fp_result.method == :FPplus
    @test fp_result.fractionality == 0.0
    @test all(feasible_integer(ctx, round.(Int, x)) for x in fp_result.solutions)
    @test length(fp_result.solutions) == length(fp_result.objectives)

    # Fractional LB vectors select FPGPR; check the returned solution/objective
    # pairing and Pareto property without asserting a stochastic solution count.
    fpgpr_result = FPGPR(ctx, [[0.25, 0.0], [1.0, 0.25]], 0.25)
    @test fpgpr_result.method == :FPGPR
    @test fpgpr_result.fractionality == 1.0
    @test length(fpgpr_result.solutions) == length(fpgpr_result.objectives)
    @test all(feasible_integer(ctx, round.(Int, x)) for x in fpgpr_result.solutions)
    @test all(isapprox(fpgpr_result.objectives[i],
                       getobjval(fpgpr_result.solutions[i], ctx))
              for i in eachindex(fpgpr_result.solutions))
    @test all(!dominated(fpgpr_result.objectives[i],
                         [fpgpr_result.objectives[j]
                          for j in eachindex(fpgpr_result.objectives) if j != i])
              for i in eachindex(fpgpr_result.objectives))
end

@testset "Degenerate path-relinking input" begin
    B = zeros(0, 1)
    ctx = build_context(B, Float64[], String[], ones(3, 1))
    X, Y, added = GPR(ctx, [[1.0]], [[1.0, 1.0, 1.0]], 0.1)
    @test length(X) == length(Y) == 1
    @test added == 0
end
