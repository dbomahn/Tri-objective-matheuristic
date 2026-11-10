# Code architecture

## Algorithm components

The reorganised source separates three concerns:

1. **Input and lower-bound data** — `Data.jl` reads the MPS/VLP model representation and stored Bensolve variable/objective data.
2. **Feasibility-pump search** — `FeasibilityPump.jl` contains the fractional feasibility-pump routine (`fractionalFP`), feasibility checks, LP projection, bit-flipping/tabu helpers and solution post-processing.
3. **Path relinking** — `GenericPathRelinking.jl` contains neighbourhood generation, selection of the next guiding solution, generic path relinking (`GPR`) and the combined `FPGPR` routine.

The algorithm files retain the legacy implementation's dependency on shared JuMP model objects and global names (`mip`, `mx`, `dist` and `x`). They are separated source files, not yet a self-contained Julia module with an independently validated API.

## FP⁺ versus FPGPR

Section 4.1 of the paper selects the method according to the fractionality of Bensolve's LB solutions. In the computational study, `allowedFractionality = 20%`:

- If measured fractionality is at most 20%, use FP⁺.
- If it exceeds 20%, use FPGPR.

FP⁺ extends the feasibility-pump search by generating additional fractional solutions. FPGPR combines a basic feasibility pump with generic path relinking. The instance family is informative about expected fractionality, but should not replace the measured runtime criterion.

The reported MOAP LB solutions are integral (the assignment formulation is totally unimodular), while MOKP and TOMIPLIB exhibit high fractionality in the reported experiments. TOFLP should be assessed using its measured LB-solution fractionality rather than assigned a method solely by family name.

## Current integration status

There is not yet a single experiment driver that performs Bensolve data loading, calculates fractionality, dispatches to FP⁺ or FPGPR, applies family-specific time limits and records all computational outputs. Historical drivers are retained under `scripts/legacy/` for audit and provenance, not as validated run commands.
