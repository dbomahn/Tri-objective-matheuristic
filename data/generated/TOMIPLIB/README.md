# TOMIPLIB benchmark instances

## Origin and selection

The TOMIPLIB set consists of newly generated tri-objective binary integer linear programming instances based on single-objective binary integer linear programming instances from [MIPLIB 2017](https://miplib.zib.de/).

The construction described in the paper is:

1. Select 90 single-objective binary integer linear programming instances categorised as *easy*, meaning that they can be solved within one hour by a single-objective solver using its default settings.
2. Generate coefficients for the second and third objective functions from the original objective coefficients, using the coefficient-generation scheme below. This is similar to the approach of Adelgren and Gupte (2022).
3. The resulting tri-objective instances are computationally difficult. For the computational study, retain the five instances for which Bensolve can find the complete set of lower-bound (LB) solutions within 15 minutes.

Thus, the computational benchmark archive contains **five instances**, selected from the initial set of 90 candidates.

## Additional objective coefficient generation

Let (c^1=(c^1_1,ldots,c^1_n)) be the original objective coefficient vector, and let (c^2) and (c^3) be the coefficient vectors for the two additional objectives.

- If (n>1), shuffle the entries of (c^1) to obtain two different coefficient vectors, (c^2) and (c^3). In other words, the added objective coefficient vectors are different permutations of the original vector.
- If (n=1), generate (c^2_1) and (c^3_1) independently according to the sign of the original coefficient:
  - If (c^1_1geq 0), draw each from the integer interval ([1,10]).
  - If (c^1_1<0), draw each from the integer interval ([-10,-1]).

The published description refers to shuffling the original coefficients twice to obtain two different vectors. It does not specify the random seed or the exact procedure used to ensure distinct permutations. Consequently, the documentation describes the stated scheme but does not claim that the original instances can be regenerated identically.

## Computational selection criterion

An instance is included in the five-instance TOMIPLIB set if Bensolve can compute its complete LB solution set within **15 minutes**. This criterion is distinct from the one-hour single-objective solve-time criterion used to select the initial 90 MIPLIB candidates.

## Files

The five instances used in the computational study are provided in [`TOMIPLIB.zip`](../TOMIPLIB.zip). The archive is retained as supplied.

## References

- MIPLIB 2017: [https://miplib.zib.de/](https://miplib.zib.de/).
- Adelgren, N. and Gupte, A. (2022). The objective-coefficient generation approach is cited in the paper as a similar method; see the paper's bibliography for the full reference.
- An, D., Parragh, S. N., Sinnl, M. and Tricoire, F. (2024). *A matheuristic for tri-objective binary integer linear programming*. Computers & Operations Research, 161, 106397. [https://doi.org/10.1016/j.cor.2023.106397](https://doi.org/10.1016/j.cor.2023.106397).
