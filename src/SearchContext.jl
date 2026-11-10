"""
Solver state shared by feasibility-pump and path-relinking routines.
The integer and relaxation models must encode the same feasible region.
C is a 3-by-n objective coefficient matrix; all three objectives are maximised.
"""
struct SearchContext{M1,M2,V1,V2,R}
    integer_model::M1
    relaxation_model::M2
    binary_variables::V1
    relaxation_variables::V2
    C::Matrix{Float64}
    rng::R
    function SearchContext(integer_model, relaxation_model, binary_variables,
                           relaxation_variables, C::AbstractMatrix; rng=Random.default_rng())
        size(C, 1) == 3 || throw(ArgumentError("Exactly three objective rows are required."))
        size(C, 2) == length(binary_variables) ||
            throw(DimensionMismatch("Objective matrix and variable dimensions differ."))
        length(binary_variables) == length(relaxation_variables) ||
            throw(DimensionMismatch("Integer and relaxation model dimensions differ."))
        new{typeof(integer_model),typeof(relaxation_model),typeof(binary_variables),
            typeof(relaxation_variables),typeof(rng)}(
            integer_model, relaxation_model, binary_variables, relaxation_variables,
            Matrix{Float64}(C), rng)
    end
end

nvariables(ctx::SearchContext) = length(ctx.binary_variables)
getobjval(x, C::AbstractMatrix) = vec(C * x)
getobjval(x, ctx::SearchContext) = getobjval(x, ctx.C)

"""Check feasibility while always removing temporary variable fixings afterwards."""
function feasible_integer(ctx::SearchContext, candidate::AbstractVector{<:Real})
    length(candidate) == nvariables(ctx) ||
        throw(DimensionMismatch("Candidate has the wrong number of variables."))
    all(v -> isapprox(v, round(v); atol=1e-8), candidate) ||
        throw(ArgumentError("Candidate must be integer-valued."))
    fixed = VariableRef[]
    try
        for (v, value_) in zip(ctx.binary_variables, candidate)
            fix(v, round(Int, value_); force=true)
            push!(fixed, v)
        end
        optimize!(ctx.integer_model)
        return has_values(ctx.integer_model) &&
            primal_status(ctx.integer_model) in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    finally
        foreach(v -> is_fixed(v) && unfix(v), fixed)
    end
end

"""LP projection of a rounded vector using its Hamming-distance objective."""
function fbsearch(ctx::SearchContext, rounded::AbstractVector{<:Integer})
    length(rounded) == nvariables(ctx) ||
        throw(DimensionMismatch("Candidate has the wrong number of variables."))
    x = ctx.relaxation_variables
    idx0, idx1 = findall(==(0), rounded), findall(==(1), rounded)
    @objective(ctx.relaxation_model, Min,
        sum(x[i] for i in idx0) + sum(1.0 - x[j] for j in idx1))
    optimize!(ctx.relaxation_model)
    return has_values(ctx.relaxation_model) ? value.(x) : nothing
end
