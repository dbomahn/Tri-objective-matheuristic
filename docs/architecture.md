# Architecture and algorithm flow

## Solver context

The search routines receive a `SearchContext` rather than relying on global names. It stores the JuMP mixed-integer feasibility model, the continuous relaxation model, both sets of variable references, the 3-by-n objective coefficient matrix and a random-number generator.

`build_context(B, rhs, senses, C)` builds both models from the same constraint system. Row senses are `"l"` for >=, `"u"` for <=, and `"s"` for equality. GLPK is the default optimizer; another JuMP-compatible optimizer can be supplied.

## Algorithm components

- `Data.jl`: model construction from a constraint matrix.
- `SearchContext.jl`: feasibility checking and LP projection.
- `FeasibilityPump.jl`: fractional feasibility-pump search, tabu perturbations and Pareto filtering.
- `GenericPathRelinking.jl`: neighbourhood generation, next-solution selection, GPR and the FPGPR two-phase search.
- `Matheuristic.jl`: per-instance FP⁺/FPGPR dispatch.

All solver models and variable references are passed explicitly. The reorganised core has no dependency on global `mip`, `mx`, `dist` or `x` bindings.

## FP⁺ versus FPGPR

Section 4.1 selects the method according to the fractionality of Bensolve's lower-bound solutions. The computational study uses `allowedFractionality = 20%`. The implementation measures the fraction of LB solution vectors containing at least one fractional variable:

- fractionality <= 20%: FP⁺;
- fractionality > 20%: FPGPR.

The instance family is informative about expected fractionality, but dispatch is based on the measured value for each instance.

## Algorithm flow and time budget

- **FP⁺** runs the fractional feasibility-pump search on the LB solutions.
- **FPGPR** performs path relinking with feasibility-pump intensification around each neighbour, followed by a pure GPR phase. Both phases share the single algorithm time limit.

The default iteration caps follow the historical implementation (10 pump iterations and 20 path-relinking iterations). The experiment driver must supply the benchmark-family-specific time limit from the paper.

## Remaining integration

The core search API and model builder are separated from file I/O. MPS/VLP/JLD2 readers and a paper-specific experiment driver still need to be connected to the actual instance and solution formats. This is separate from the algorithm-level defects fixed in the new core.
