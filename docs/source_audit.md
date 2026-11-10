# Source audit

## Sources reviewed

The reorganisation is based on these files in the private `dbomahn/phd` repository:

- `codes/FPGPRfunctions.jl` — principal shared function library.
- `codes/(F)FPGPR.jl` — standalone experiment script containing duplicated functions, model construction and experiment execution.
- `codes/FP+GPR.jl` — sequential FP/GPR experiment script.
- `codes/kpFPGPR.jl` — separate knapsack-specific implementation.
- `codes/InstanceGenerator.jl` — unrelated/uncorroborated instance-processing script.

The shared function library was separated into `src/Data.jl`, `src/FeasibilityPump.jl` and `src/GenericPathRelinking.jl`. Historical scripts are retained under `scripts/legacy/`.

## Known issues requiring reconciliation

These are source-inspection findings, not runtime test results:

1. The original `FPGPRfunctions.jl` contains top-level `ARGS[1:5]` accesses, unsuitable for a reusable library. These were removed from the reorganised data/functions files.
2. The library's `GPR` routine uses `iter` in its loop condition and increments it without initialising it in that function.
3. `flipoper` contains a reference to `r` in the multi-bit-flip loop where the loop variable is `i`.
4. Search routines depend on globally defined JuMP models and variables (`mip`, `mx`, `dist`, `x`); the current split does not encapsulate these dependencies.
5. The standalone `(F)FPGPR.jl` script references `fpY` although the preceding `fractionalFP` call returns only `fpX` and a count. It also relies on positional command-line arguments and hard-coded output paths.
6. `FP+GPR.jl` includes `FFPGPRfunctions.jl`, absent from the inspected tree (the available filename is `FPGPRfunctions.jl`). It calls `weightFP`, not found in the inspected shared library; its GPR calls use an argument order inconsistent with the function signature; and it refers to `fpY` although it is not produced by the preceding fractional-pump call.
7. The knapsack-specific implementation has its own data/model assumptions and should not be merged into the generic MPS/VLP implementation without a separate compatibility review.

## InstanceGenerator.jl

The inspected file points to a hard-coded Windows directory named `tri-objBoland2015b`. Its only generation loop is commented out and would append random third-objective coefficients to files in that directory. The active code reads one file from the second directory entry, formats objective/constraint data and constructs a CPLEX model; it does not generate or save a TOFLP or TOMIPLIB instance.

This does not match the published TOFLP construction (Gadegaard et al. single-source capacitated facility-location instances plus demand-coverage objective) or the published TOMIPLIB construction (MIPLIB 2017 instances with shuffled objective coefficients). The file is therefore not presented as the generator for either benchmark family.

## Validation status

No Julia execution, solver run or numerical comparison was performed during this source reorganisation. The legacy drivers are retained as historical source, not certified as executable.
