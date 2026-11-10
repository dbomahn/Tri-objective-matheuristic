function flip(x::AbstractVector{<:Integer}, position::Int)
    y = copy(x)
    y[position] = 1 - y[position]
    return y
end

"""Find an unvisited one-bit or multi-bit perturbation, or return nothing."""
function flipoper(ctx::SearchContext, tabu, fractional, rounded)
    n = length(rounded)
    n == 0 && return nothing
    for j in sortperm(abs.(fractional .- rounded); rev=true)
        candidate = flip(rounded, j)
        candidate in tabu || return candidate
    end
    if n > 1
        for _ in 1:min(n, 10)
            count = rand(ctx.rng, ceil(Int, n / 2):(n - 1))
            positions = randperm(ctx.rng, n)[1:count]
            candidate = copy(rounded)
            for j in positions
                candidate[j] = 1 - candidate[j]
            end
            candidate in tabu || return candidate
        end
    end
    return nothing
end

"""True if x is weakly dominated by a point in P (maximisation)."""
dominated(x, P) = any(p -> all(x .<= p), P)

"""Pareto-filter maximisation points while preserving solution/objective pairing."""
function pareto_filter(solutions, objectives)
    length(solutions) == length(objectives) ||
        throw(DimensionMismatch("Solution and objective collections must have equal length."))
    keep = trues(length(objectives))
    for i in eachindex(objectives)
        keep[i] || continue
        for j in eachindex(objectives)
            (i == j || !keep[j]) && continue
            oi, oj = objectives[i], objectives[j]
            if all(oi .>= oj) && any(oi .> oj)
                keep[j] = false
            elseif all(oj .>= oi) && any(oj .> oi)
                keep[i] = false
                break
            elseif oi == oj && j > i
                keep[j] = false
            end
        end
    end
    idx = findall(keep)
    return solutions[idx], objectives[idx]
end

"""
Fractional feasibility-pump search from the Bensolve LB decision vectors.
Returns all starting and newly found feasible integer vectors, plus the number
of newly found vectors. Their objective values are recomputed by the caller.
"""
function fractionalFP(ctx::SearchContext, initial_solutions, time_limit::Real;
                      max_iterations::Int=10)
    time_limit >= 0 || throw(ArgumentError("time_limit must be non-negative."))
    X = [Float64.(x) for x in initial_solutions]
    length(X) < 2 && return X, 0
    tabu = Vector{Vector{Int}}()
    explored_pairs = Set{Tuple{Int,Int}}()
    new_solutions = 0
    deadline = time() + time_limit

    while time() < deadline && length(X) >= 2
        i, j = sort(randperm(ctx.rng, length(X))[1:2])
        pair = (i, j)
        if pair in explored_pairs
            length(explored_pairs) >= binomial(length(X), 2) && break
            continue
        end
        push!(explored_pairs, pair)
        lambda = rand(ctx.rng, 0:10) / 10
        target = lambda .* X[i] .+ (1 - lambda) .* X[j]
        for _ in 1:max_iterations
            time() >= deadline && break
            rounded = round.(Int, target)
            already_present = any(x -> round.(Int, x) == rounded, X)
            if !already_present && feasible_integer(ctx, rounded)
                push!(X, Float64.(rounded))
                new_solutions += 1
                break
            end
            if rounded in tabu
                alternative = flipoper(ctx, tabu, target, rounded)
                alternative === nothing && break
                rounded = alternative
                already_present = any(x -> round.(Int, x) == rounded, X)
                if !already_present && feasible_integer(ctx, rounded)
                    push!(X, Float64.(rounded))
                    new_solutions += 1
                    break
                end
            end
            push!(tabu, copy(rounded))
            target = fbsearch(ctx, rounded)
            target === nothing && break
        end
    end
    return X, new_solutions
end

FPplus(ctx::SearchContext, initial_solutions, time_limit::Real; kwargs...) =
    fractionalFP(ctx, initial_solutions, time_limit; kwargs...)
