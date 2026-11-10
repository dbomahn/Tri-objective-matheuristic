"""
Compute the fraction of Bensolve LB solution vectors that contain at least one
fractional variable. This is the per-instance measure used for dispatch.
"""
function select_algorithm(lb_solutions; allowed_fractionality::Real=0.20)
    0 <= allowed_fractionality <= 1 ||
        throw(ArgumentError("allowed_fractionality must lie in [0,1]."))
    isempty(lb_solutions) && throw(ArgumentError("At least one LB solution is required."))
    n = length(first(lb_solutions))
    n > 0 || throw(ArgumentError("LB solutions must contain variables."))
    all(length(x) == n for x in lb_solutions) ||
        throw(DimensionMismatch("LB solutions have inconsistent dimensions."))
    fractional_solutions = count(x -> any(v -> 1e-6 < v < 1 - 1e-6, x), lb_solutions)
    fractionality = fractional_solutions / length(lb_solutions)
    return fractionality <= allowed_fractionality ? :FPplus : :FPGPR, fractionality
end

"""
Run the fractionality-based dispatch described in Section 4.1.
For FP⁺, run the fractional pump for the full algorithm budget. For FPGPR,
split the budget between the pump and generic path relinking.
"""
function FPGPR(ctx::SearchContext, lb_solutions, time_limit::Real;
               allowed_fractionality::Real=0.20, pump_fraction::Real=0.5)
    time_limit >= 0 || throw(ArgumentError("time_limit must be non-negative."))
    0 <= pump_fraction <= 1 || throw(ArgumentError("pump_fraction must lie in [0,1]."))
    method, fractionality = select_algorithm(lb_solutions;
        allowed_fractionality=allowed_fractionality)
    pump_budget = method == :FPplus ? time_limit : time_limit * pump_fraction
    X, added_by_pump = FPplus(ctx, lb_solutions, pump_budget)
    Y = [getobjval(x, ctx) for x in X]
    X, Y = pareto_filter(X, Y)

    if method == :FPplus || time_limit <= pump_budget
        return (method=method, fractionality=fractionality, solutions=X,
                objectives=Y, pump_added=added_by_pump, relinking_added=0)
    end

    X, Y, added_by_relinking = GPR(ctx, X, Y, time_limit - pump_budget)
    return (method=method, fractionality=fractionality, solutions=X,
            objectives=Y, pump_added=added_by_pump,
            relinking_added=added_by_relinking)
end
