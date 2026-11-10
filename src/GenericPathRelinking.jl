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

"""
Generic path relinking. The iteration counter is local to each path and the
routine safely returns when fewer than two starting solutions are available.
"""
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
        pairs = [(i, j) for i in 1:length(X)-1 for j in i+1:length(X)
                 if !((i, j) in tried_pairs)]
        isempty(pairs) && break
        pair = rand(ctx.rng, pairs)
        push!(tried_pairs, pair)
        i, j = pair
        current = copy(X[i])
        guiding = copy(X[j])

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
