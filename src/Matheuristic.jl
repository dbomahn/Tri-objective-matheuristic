"""
Compute the fraction of Bensolve LB solution vectors containing at least one
fractional variable. This is the per-instance measure used for method dispatch.
"""
function select_algorithm(lb_solutions; allowed_fractionality::Real=0.20)
    0 <= allowed_fractionality <= 1 ||
        throw(ArgumentError("allowed_fractionality must lie in [0,1]."))
    isempty(lb_solutions) && throw(ArgumentError("At least one LB solution is required."))
    n = length(first(lb_solutions))
    n > 0 || throw(ArgumentError("LB solutions must contain variables."))
    all(length(x) == n for x in lb_solutions) ||
        throw(DimensionMismatch("LB solutions have inconsistent dimensions."))
    fractional_solutions = count(
        x -> any(v -> 1e-6 < v < 1 - 1e-6, x), lb_solutions)
    fractionality = fractional_solutions / length(lb_solutions)
    return fractionality <= allowed_fractionality ? :FPplus : :FPGPR, fractionality
end

"""
Dispatch according to the measured LB-solution fractionality and run the
selected algorithm under the supplied algorithm time limit.
"""
function FPGPR(ctx::SearchContext, lb_solutions, time_limit::Real;
               allowed_fractionality::Real=0.20)
    time_limit >= 0 || throw(ArgumentError("time_limit must be non-negative."))
    method, fractionality = select_algorithm(lb_solutions;
        allowed_fractionality=allowed_fractionality)

    if method == :FPplus
        X, added = FPplus(ctx, lb_solutions, time_limit)
        Y = [getobjval(x, ctx) for x in X]
        X, Y = pareto_filter(X, Y)
        return (method=method, fractionality=fractionality, solutions=X,
                objectives=Y, new_solutions=added)
    end

    initial_objectives = [getobjval(x, ctx) for x in lb_solutions]
    X, Y, added = FPGPR_search(ctx, lb_solutions, initial_objectives, time_limit)
    return (method=method, fractionality=fractionality, solutions=X,
            objectives=Y, new_solutions=added)
end
