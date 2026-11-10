module TriObjectiveMatheuristic
using JuMP, LinearAlgebra, Random, StatsBase, GLPK
import MathOptInterface as MOI
include("Data.jl")
include("SearchContext.jl")
include("FeasibilityPump.jl")
include("GenericPathRelinking.jl")
include("Matheuristic.jl")
export SearchContext, build_context, feasible_integer, fbsearch, fractionalFP,
       FPplus, GPR, FPGPR, getobjval, pareto_filter, select_algorithm, dominated
end
