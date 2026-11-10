"""Generate unvisited neighbours along the differing coordinates."""
function createNB(ctx::SearchContext, solution, differing_indices, explored)
    neighbours = Vector{Vector{Float64}}()
    objectives = Vector{Vector{Float64}}()
    for i in differing_indices
        values = isapprox(solution[i], round(solution[i]); atol=1e-8) ?
            (1.0 - round(solution[i]),) : (1.0, 0.0)
        for value_ in values
            candidate = Float64.(solution)
            candidate[i] = value_
            if !(candidate in explored)
                push!(neighbours, candidate)
                push!(objectives, getobjval(candidate, ctx))
            end
        end
    end
    return neighbours, objectives
end

"""Select the next solution by summed ranks of its three objective ratios."""
function nextSI(ctx::SearchContext, neighbours, neighbour_objectives, current)
    isempty(neighbours) && return nothing
    length(neighbours) == 1 && return only(neighbours)
    current_objectives = getobjval(current, ctx)
    ratios = Matrix{Float64}(undef, length(neighbour_objectives), 3)
    for i in eachindex(neighbour_objectives), k in 1:3
        denominator = current_objectives[k]
        ratio = abs(denominator) > eps(Float64) ?
            neighbour_objectives[i][k] / denominator :
            neighbour_objectives[i][k] - current_objectives[k]
        ratios[i, k] = isfinite(ratio) ? ratio :
            neighbour_objectives[i][k] - current_objectives[k]
    end
    ranks = hcat((tiedrank(ratios[:, k]) for k in 1:3)...)
    totals = vec(sum(ranks; dims=2))
    best = findall(==(maximum(totals)), totals)
    return neighbours[rand(ctx.rng, best)]
end

function next_untried_pair(ctx::SearchContext, n::Int, tried_pairs)
    n < 2 && return nothing
    available = [(i, j) for i in 1:(n - 1) for j in (i + 1):n
                 if !((i, j) in tried_pairs)]
    isempty(available) && return nothing
    return rand(ctx.rng, available)
end

"""Generic path relinking without the feasibility-pump intensification phase."""
function GPR(ctx::SearchContext, initial_solutions, initial_objectives, time_limit::Real;
             max_path_iterations::Int=20)
    time_limit >= 0 || throw(ArgumentError("time_limit must be non-negative."))
    length(initial_solutions) == length(initial_objectives) ||
        throw(DimensionMismatch("Solution and objective collections must have equal length."))
    X = [Float64.(x) for x in initial_solutions]
    Y = [Float64.(y) for y in initial_objectives]
    length(X) < 2 && return X, Y, 0

    deadline = time() + time_limit
    explored = Vector{Vector{Float64}}()
    tried_pairs = Set{Tuple{Int,Int}}()
    added = 0

    while time() < deadline && length(X) >= 2
        pair = next_untried_pair(ctx, length(X), tried_pairs)
        pair === nothing && break
        push!(tried_pairs, pair)
        i, j = pair
        current, guiding = copy(X[i]), copy(X[j])

        for _ in 1:max_path_iterations
            time() >= deadline && break
            all(k -> isapprox(current[k], guiding[k]; atol=1e-8), eachindex(current)) && break
            differing = findall(k -> !isapprox(current[k], guiding[k]; atol=1e-8),
                                eachindex(current))
            neighbours, neighbour_objectives = createNB(ctx, current, differing, explored)
            isempty(neighbours) && break

            for k in eachindex(neighbours)
                time() >= deadline && break
                candidate = round.(Int, neighbours[k])
                present = any(x -> round.(Int, x) == candidate, X)
                if !present && feasible_integer(ctx, candidate)
                    push!(X, Float64.(candidate))
                    push!(Y, getobjval(candidate, ctx))
                    added += 1
                end
            end

            next = nextSI(ctx, neighbours, neighbour_objectives, current)
            next === nothing && break
            push!(explored, copy(next))
            current = next
        end
    end
    X, Y = pareto_filter(X, Y)
    return X, Y, added
end

"""
FPGPR implementation: path relinking with feasibility-pump intensification
around each neighbour, followed by a pure generic path-relinking phase.
Both phases share the supplied time limit.
"""
function FPGPR_search(ctx::SearchContext, initial_solutions, initial_objectives,
                      time_limit::Real; gpr_iterations::Int=20, fp_iterations::Int=10)
    time_limit >= 0 || throw(ArgumentError("time_limit must be non-negative."))
    X = [Float64.(x) for x in initial_solutions]
    Y = [Float64.(y) for y in initial_objectives]
    length(X) == length(Y) || throw(DimensionMismatch("Solution/objective lengths differ."))
    length(X) < 2 && return X, Y, 0

    deadline = time() + time_limit
    explored = Vector{Vector{Float64}}()
    tabu = Vector{Vector{Int}}()
    tried_pairs = Set{Tuple{Int,Int}}()
    added = 0
    initial_pool_size = length(X)

    # Phase 1: intensify paths between the original Bensolve LB solutions.
    while time() < deadline
        pair = next_untried_pair(ctx, initial_pool_size, tried_pairs)
        pair === nothing && break
        push!(tried_pairs, pair)
        i, j = pair
        current, guiding = copy(X[i]), copy(X[j])
        path_iteration = 0

        while any(k -> !isapprox(current[k], guiding[k]; atol=1e-8), eachindex(current)) &&
              path_iteration < gpr_iterations && time() < deadline
            differing = findall(k -> !isapprox(current[k], guiding[k]; atol=1e-8),
                                eachindex(current))
            neighbours, neighbour_objectives = createNB(ctx, current, differing, explored)
            isempty(neighbours) && break

            for l in eachindex(neighbours)
                target = copy(neighbours[l])
                search_done = false
                for _ in 1:fp_iterations
                    time() >= deadline && break
                    rounded = round.(Int, target)
                    objective = getobjval(rounded, ctx)
                    present = any(y -> y == objective, Y)

                    if !present && feasible_integer(ctx, rounded)
                        push!(X, Float64.(rounded))
                        push!(Y, objective)
                        added += 1
                        search_done = true
                        break
                    end

                    if rounded in tabu
                        alternative = flipoper(ctx, tabu, target, rounded)
                        if alternative === nothing
                            search_done = true
                            break
                        end
                        rounded = alternative
                        objective = getobjval(rounded, ctx)
                        present = any(y -> y == objective, Y)
                        if !present && feasible_integer(ctx, rounded)
                            push!(X, Float64.(rounded))
                            push!(Y, objective)
                            added += 1
                            search_done = true
                            break
                        end
                    end

                    if !search_done
                        push!(tabu, copy(rounded))
                        projected = fbsearch(ctx, rounded)
                        if projected === nothing
                            search_done = true
                            break
                        end
                        target = projected
                    end
                end
            end

            next = nextSI(ctx, neighbours, neighbour_objectives, current)
            next === nothing && break
            push!(explored, copy(next))
            current = next
            path_iteration += 1
        end
    end

    # Phase 2: pure GPR on the expanded candidate set.
    while time() < deadline
        pair = next_untried_pair(ctx, length(X), tried_pairs)
        pair === nothing && break
        push!(tried_pairs, pair)
        i, j = pair
        current, guiding = copy(X[i]), copy(X[j])
        path_iteration = 0

        while any(k -> !isapprox(current[k], guiding[k]; atol=1e-8), eachindex(current)) &&
              path_iteration < gpr_iterations && time() < deadline
            differing = findall(k -> !isapprox(current[k], guiding[k]; atol=1e-8),
                                eachindex(current))
            neighbours, neighbour_objectives = createNB(ctx, current, differing, explored)
            isempty(neighbours) && break

            for l in eachindex(neighbours)
                time() >= deadline && break
                candidate = round.(Int, neighbours[l])
                present = any(x -> round.(Int, x) == candidate, X)
                if !present && feasible_integer(ctx, candidate)
                    push!(X, Float64.(candidate))
                    push!(Y, getobjval(candidate, ctx))
                    added += 1
                end
            end

            next = nextSI(ctx, neighbours, neighbour_objectives, current)
            next === nothing && break
            push!(explored, copy(next))
            current = next
            path_iteration += 1
        end
    end

    X, Y = pareto_filter(X, Y)
    return X, Y, added
end
