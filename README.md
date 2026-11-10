# Tri-objective matheuristic

Julia implementation accompanying **A matheuristic for tri-objective binary integer linear programming**.

This repository separates algorithm source from historical experiment scripts and documents the benchmark data used in the computational study. It is a source-level reorganisation of code from the author's PhD repository; it is not yet a verified, turnkey reproduction package.

## Repository layout

- `src/Data.jl`: input-data structures and readers.
- `src/FeasibilityPump.jl`: feasibility-pump routines and shared solution-processing helpers.
- `src/GenericPathRelinking.jl`: generic path-relinking routines and the combined FPGPR routine.
- `scripts/legacy/`: original experiment scripts retained for provenance; these are not advertised as working entry points.
- `data/`: locations for the benchmark instances, which will be added separately.
- `docs/`: architecture, benchmark provenance and source audit.

All Julia source files use the `.jl` extension.

## Algorithm selection

The paper's FP⁺/FPGPR choice is based on the **fractionality of the lower-bound (LB) solutions returned by Bensolve**, not simply on a fixed instance-family label. The computational study uses an allowed-fractionality threshold of 20%: at or below the threshold, use FP⁺; above it, use FPGPR. See [docs/architecture.md](docs/architecture.md).

## Benchmark data

Instance files are not included in this initial reorganisation. Their intended locations and published construction methodology are described in [docs/instance_generation.md](docs/instance_generation.md). The original generation scripts are unavailable, so this repository does not claim to regenerate those files.

## Status and limitations

The source has been separated into data, feasibility-pump and path-relinking files. Historical drivers are retained separately and marked as legacy. Several source-level issues remain unresolved; see [docs/source_audit.md](docs/source_audit.md).

No computational experiment or Julia test suite has been run. A reproducible experiment driver, dependency manifest and automated tests should be added after the historical implementation has been reconciled and the benchmark files are uploaded.

## Paper

An, D., Lee, C. and Parragh, S. N. (2024). *A matheuristic for tri-objective binary integer linear programming*. Computers & Operations Research, 161, 106397. https://doi.org/10.1016/j.cor.2023.106397
