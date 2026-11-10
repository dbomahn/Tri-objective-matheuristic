"""
Build the two JuMP models required by the search from a constraint matrix.
Each row sense is one of "l" (>=), "u" (<=), or "s" (=).
This keeps solver models explicit and avoids global model/variable references.
"""
function build_context(B::AbstractMatrix, rhs::AbstractVector, senses::AbstractVector,
                       C::AbstractMatrix; optimizer=GLPK.Optimizer,
                       rng=Random.default_rng(), silent::Bool=true)
    m, n = size(B)
    length(rhs) == m || throw(DimensionMismatch("RHS length must equal the row count of B."))
    length(senses) == m || throw(DimensionMismatch("Sense count must equal the row count of B."))
    size(C) == (3, n) || throw(DimensionMismatch("C must be a 3-by-n objective matrix."))
    all(s -> s in ("l", "u", "s"), senses) ||
        throw(ArgumentError("Constraint senses must be l, u, or s."))

    integer_model = Model(optimizer)
    relaxation_model = Model(optimizer)
    silent && (set_silent(integer_model); set_silent(relaxation_model))
    @variable(integer_model, z[1:n], Bin)
    @variable(relaxation_model, 0 <= x[1:n] <= 1)

    for i in 1:m
        lhs_integer = sum(B[i, j] * z[j] for j in 1:n)
        lhs_relaxation = sum(B[i, j] * x[j] for j in 1:n)
        if senses[i] == "l"
            @constraint(integer_model, lhs_integer >= rhs[i])
            @constraint(relaxation_model, lhs_relaxation >= rhs[i])
        elseif senses[i] == "u"
            @constraint(integer_model, lhs_integer <= rhs[i])
            @constraint(relaxation_model, lhs_relaxation <= rhs[i])
        else
            @constraint(integer_model, lhs_integer == rhs[i])
            @constraint(relaxation_model, lhs_relaxation == rhs[i])
        end
    end
    @objective(integer_model, Min, 0)
    @objective(relaxation_model, Min, 0)
    return SearchContext(integer_model, relaxation_model, z, x, C; rng=rng)
end
