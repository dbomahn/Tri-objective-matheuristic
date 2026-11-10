# Architecture and algorithm flow

## Solver context

The search routines receive a `SearchContext` rather than relying on global names. It stores:
- the JuMP mixed-integer feasibility model;
- the JuMP continuous relaxation model;
- the corresponding integer and relaxed variable references;
- the 3-by-n objective coefficient matrix; and
- a random-number generator.

`build_context(B, rhs, senses, C)` builds both models from the same constraint system. Row senses are `"l"` for >=, `"u"` for <=, and `"s"` for equality. GLPK is the default optimizer; another JuMP-compatible optimizer can be supplied.

## Algorithm components

- `Data.jl`: model construction from a constraint matrix.
- `SearchContext.jl`: feasibility checking and LP projection.
- `FeasibilityPump.jl`: fractional feasibility-pump search, tabu perturbations and Pareto filtering.
- `GenericPathRelinking.jl`: neighbourhood generation, next-solution selection and GPR.
- `Matheuristic.jl`: per-instance FP⁺/FPGPR dispatch.

All solver models and variable references are passed explicitly. There are no dependencies on global `mip`, `mx`, `dist` or `x` bindings in the reorganised core.

## FP⁺ versus FPGPR

Section 4.1 selects the method according to the fractionality of Bensolve's lower-bound solutions. The computational study uses `allowedFractionality = 20%`. The implementation measures the fraction of LB solution vectors containing at least one fractional variable:

- fractionality <= 20%: FP⁺;
- fractionality > 20%: FPGPR.

The instance family is informative about expected fractionality, but the dispatch is based on the measured value for the instance.

## Time budgets

The `FPGPR` dispatcher accepts the algorithm time budget from the experiment protocol. For FP⁺, the full budget is assigned to the pump. For FPGPR, the default split is 50% for the pump and the remainder for path relinking; this is configurable through `pump_fraction`. Benchmark-specific budgets should be passed by the experiment driver in accordance with the paper.

## Remaining integration

The core search API and model builder are separated from file I/O. The MPS/VLP/JLD2 readers and paper-specific experiment driver still need to be connected to the actual instance/solution formats. This is distinct from the fixed algorithm-level defects in the old scripts.
