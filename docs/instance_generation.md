# Benchmark instance provenance

The repository includes the TOFLP and TOMIPLIB files used in the computational study. The original instance-generation scripts are not available. The descriptions below document the published construction process; they are not a claim that this repository can regenerate byte-identical instances.

## MOAP

The multi-objective assignment problem instances are sourced from Kirlik and Sayın (2014). The paper uses 10 subclasses with 10 instances each (100 instances). These are externally hosted and are not copied into this repository. See [data/external/README.md](../data/external/README.md) for the publication and benchmark-library links.

## MOKP

The multi-objective knapsack problem instances are sourced from Kirlik and Sayın (2014). The paper uses 10 subclasses with 10 instances each (100 instances). These are externally hosted and are not copied into this repository. See [data/external/README.md](../data/external/README.md) for the publication and benchmark-library links.

## TOFLP

The tri-objective facility-location instances extend the bi-objective single-source capacitated facility-location model of Gadegaard et al. (2019) by adding total covered customer demand as the third objective. The model, assumptions, instance naming convention and file layout are documented in [data/generated/TOFLP/README.md](../data/generated/TOFLP/README.md).

The set contains 12 subclasses, with 10 instances per subclass (120 instances). The number of facilities ranges from 5 to 60 in increments of five; each instance has twice as many customers as facilities. The supplied archive is [data/generated/TOFLP.zip](../data/generated/TOFLP.zip).

## TOMIPLIB

The tri-objective binary integer linear programming instances are based on single-objective binary integer linear programming instances from [MIPLIB 2017](https://miplib.zib.de/).

The published construction first selects 90 instances categorised as *easy*, meaning that a single-objective solver can solve them within one hour using default settings. The second and third objective coefficient vectors are then generated from the original objective vector: for more than one variable, two different shuffles (permutations) of the original coefficients are used; for a one-variable instance, two same-sign integer coefficients are drawn from [1,10] or [-10,-1], according to the sign of the original coefficient. This approach is described as similar to Adelgren and Gupte (2022).

The tri-objective instances proved difficult. The five instances used in the computational study are those for which Bensolve found the complete set of lower-bound solutions within 15 minutes. The detailed description is in [data/generated/TOMIPLIB/README.md](../data/generated/TOMIPLIB/README.md), and the supplied archive is [data/generated/TOMIPLIB.zip](../data/generated/TOMIPLIB.zip).

## Reproducibility note

The published descriptions establish benchmark provenance and high-level construction methodology. They do not provide the exact generation scripts, random seeds, or all implementation details needed to regenerate byte-identical files. The supplied TOFLP and TOMIPLIB files should therefore be treated as the authoritative computational-study data.
