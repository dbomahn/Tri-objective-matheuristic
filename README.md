# Tri-objective matheuristic

Julia implementation accompanying **A matheuristic for tri-objective binary integer linear programming**.

The code is organised as a reusable search core. The search does not depend on top-level global JuMP models or global variable bindings.

## Structure

- `src/Data.jl`: constructs integer and LP-relaxation models from a constraint matrix, right-hand sides and row senses.
- `src/SearchContext.jl`: explicit solver state, feasibility checks and LP projection.
- `src/FeasibilityPump.jl`: fractional feasibility-pump search and solution filtering.
- `src/GenericPathRelinking.jl`: generic path relinking.
- `src/Matheuristic.jl`: fractionality-based FP⁺/FPGPR dispatch.
- `scripts/legacy/`: historical experiment drivers retained for provenance.
- `data/`: benchmark archives and references to externally hosted benchmark sets.
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

- **TOFLP:** 120 generated instances, provided in [`data/generated/TOFLP.zip`](data/generated/TOFLP.zip). See the [TOFLP model and data README](data/generated/TOFLP/README.md).
- **TOMIPLIB:** five generated instances, provided in [`data/generated/TOMIPLIB.zip`](data/generated/TOMIPLIB.zip).
- **MOAP and MOKP:** externally hosted benchmark instances from Kirlik and Sayın (2014); this repository does not duplicate those files. See [`data/external/README.md`](data/external/README.md) for the source and publication links.

The benchmark provenance and published construction methodology are described in [`docs/instance_generation.md`](docs/instance_generation.md). The original TOFLP/TOMIPLIB generation scripts are unavailable; the repository documents their published provenance rather than claiming to regenerate the files.

## Validation and remaining integration

The package test suite passes in [GitHub Actions](https://github.com/dbomahn/Tri-objective-matheuristic/actions/runs/42889138992). These tests validate the current core implementation; they do not constitute a reproduction of the paper's computational experiments.

The MPS/VLP/Bensolve file readers and a complete experiment driver still need to be integrated and validated against the supplied benchmark formats. A full benchmark run has not yet been claimed.

## Paper

An, D., Lee, C. and Parragh, S. N. (2024). *A matheuristic for tri-objective binary integer linear programming*. Computers & Operations Research, 161, 106397. https://doi.org/10.1016/j.cor.2023.106397
