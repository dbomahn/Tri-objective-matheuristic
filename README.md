# Tri-objective matheuristic

Julia implementation accompanying **A matheuristic for tri-objective binary integer linear programming**.

The code is being reorganised into a reusable solver-independent search core. The search no longer depends on top-level global JuMP models or global variable bindings.

## Structure

- `src/Data.jl`: constructs integer and LP-relaxation models from a constraint matrix, right-hand sides and row senses.
- `src/SearchContext.jl`: explicit solver state, feasibility checks and LP projection.
- `src/FeasibilityPump.jl`: fractional feasibility-pump search and solution filtering.
- `src/GenericPathRelinking.jl`: generic path relinking.
- `src/Matheuristic.jl`: fractionality-based FP⁺/FPGPR dispatch.
- `scripts/legacy/`: unmodified historical experiment drivers retained for provenance.
- `data/`: benchmark instance locations; data files will be uploaded separately.
- `docs/`: architecture, benchmark provenance and source audit.

All Julia source files use the `.jl` extension.

## Building a search context

The core accepts a constraint matrix and associated data, so the integer and relaxed models are created in one place:

```julia
using TriObjectiveMatheuristic

ctx = build_context(B, rhs, senses, C)
result = FPGPR(ctx, bensolve_lb_solutions, 120.0)
```

Here, `B` is the constraint matrix, `rhs` is its right-hand-side vector, `senses` contains `"l"` (greater-than-or-equal), `"u"` (less-than-or-equal), or `"s"` (equality), and `C` is the 3-by-n objective matrix. The example uses a 120-second algorithm budget; benchmark-specific limits must follow the paper's experimental protocol.

## Algorithm selection

The method is selected from the measured fraction of Bensolve LB solution vectors containing at least one fractional variable. The computational study uses a threshold of 20%: at or below the threshold, FP⁺ is selected; above it, FPGPR is selected. This is a per-instance measurement, not a fixed mapping from benchmark family to algorithm.

## Benchmark data

The instance files are not included yet. Their intended locations and published construction methodology are described in [docs/instance_generation.md](docs/instance_generation.md). The original TOFLP/TOMIPLIB generation scripts are unavailable; the repository documents their published provenance rather than claiming to regenerate the files.

## Validation status

The core has been refactored to remove global model references, initialise path iteration state, repair the undefined flip index, handle insufficient candidate sets, restore temporary variable fixings, and preserve solution/objective pairing during Pareto filtering.

A Julia runtime and solver execution environment were not available during this correction, so the test suite has been added but has **not yet been executed**. The MPS/VLP/Bensolve file readers and a full experiment driver remain to be integrated once the actual instance and solution file formats are available.

## Paper

An, D., Lee, C. and Parragh, S. N. (2024). *A matheuristic for tri-objective binary integer linear programming*. Computers & Operations Research, 161, 106397. https://doi.org/10.1016/j.cor.2023.106397
