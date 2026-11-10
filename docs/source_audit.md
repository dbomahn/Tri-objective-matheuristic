# Source audit and correction record

## Historical sources reviewed

The reorganisation is based on the private `dbomahn/phd` repository:

- `codes/FPGPRfunctions.jl` — principal shared function library.
- `codes/(F)FPGPR.jl` — standalone experiment script with duplicated functions, model construction and execution.
- `codes/FP+GPR.jl` — sequential FP/GPR experiment script.
- `codes/kpFPGPR.jl` — separate knapsack-specific implementation.
- `codes/InstanceGenerator.jl` — historical instance-processing script.

The original experiment scripts are retained under `scripts/legacy/` for provenance. The new `src/` files are the corrected reusable core; the legacy scripts themselves have not been rewritten.

## Defects corrected in the reusable core

1. **Global solver state:** replaced references to global `mip`, `mx`, `dist` and `x` with an explicit `SearchContext` passed to every routine.
2. **Temporary variable fixings:** feasibility checks now remove their fixings in a `finally` block, including when optimisation throws an exception.
3. **Undefined flip index:** the multi-bit perturbation now flips the sampled coordinate `j), not the undefined variable `r`.
4. **Uninitialised GPR iteration counter:** path iteration is local to each relinking path.
5. **Insufficient candidate sets:** fractional pump and GPR return safely when fewer than two starting solutions are available.
6. **Neighbourhood alignment:** explored neighbours are filtered together with their corresponding objective vectors.
7. **Pareto filtering:** solution/objective pairing is retained while dominated points and duplicate objective vectors are removed.
8. **Fractionality dispatch:** the dispatcher measures the proportion of LB solution vectors containing fractional variables and applies the paper's 20% threshold.

## Historical driver defects retained for reference

The old `(F)FPGPR.jl` and `FP+GPR.jl` scripts are not production entry points. Source inspection found, among other things, references to `fpY` without a corresponding return value, a missing include filename, inconsistent GPR argument order, hard-coded output paths, and reliance on positional command-line arguments. They are preserved unmodified so that the original experimental code remains traceable.

The separate knapsack implementation has different data/model assumptions and has not been merged into the generic core.

## InstanceGenerator.jl

The inspected file points to a hard-coded Windows directory named `tri-objBoland2015b`. Its generation loop is commented out and would append random third-objective coefficients to files in that directory. The active code reads one file from the second directory entry, formats objective/constraint data and constructs a CPLEX model; it does not generate or save a TOFLP or TOMIPLIB instance.

This does not match the published TOFLP construction (Gadegaard et al. single-source capacitated facility-location instances plus demand-coverage objective) or TOMIPLIB construction (MIPLIB 2017 instances with shuffled objective coefficients). The file is therefore not presented as the generator for either benchmark family.

## Validation status

A focused Julia test suite has been added for Pareto filtering, fractionality dispatch, feasibility fixing/unfixing, LP projection and degenerate path-relinking inputs. It has **not been executed** because no Julia runtime is available in the current environment. Consequently, the code changes are source-level corrections, not a claim of runtime validation.
