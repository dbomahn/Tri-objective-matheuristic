# Tri-objective facility location problem (TOFLP)

## Model documentation

**Model title:** *Tri-objective facility location problem model*

**Authors:**
- Duleabom An
- Sophie N. Parragh
- Markus Sinnl
- Fabien Tricoire

Affiliations, as stated in the original model documentation:
- Institute of Production and Logistics Management, Johannes Kepler University Linz, Austria
- Institute for Transport and Logistics Management, Vienna University of Economics and Business, Austria

## Mathematical formulation

The TOFLP extends the bi-objective single-source capacitated facility-location model of Gadegaard et al. (2019) with a third objective that maximises the total demand of served customers.

Let:
- (I) denote the set of candidate facilities;
- (J) denote the set of customers;
- (x_{ij}=1) if customer (j) is served by facility (i), and 0 otherwise;
- (y_i=1) if facility (i) is opened, and 0 otherwise;
- (z_j=1) if customer (j) is served, and 0 otherwise.

The three objectives, expressed as a minimisation problem, are

[
min left(
sum_{iin I}sum_{jin J} c_{ij}x_{ij},
quad
sum_{iin I} f_i y_i,
quad
-sum_{jin J} d_j z_j
ight).
]

Subject to

[
sum_{iin I}x_{ij}=z_j
qquad orall jin J,
]

[
x_{ij}leq y_i
qquad orall iin I, jin J,
]

[
x_{ij},y_i,z_jin{0,1}
qquad orall iin I, jin J.
]

The third objective is the maximisation of covered customer demand, represented as minimisation of its negative.

## Parameters and assumptions

The original model documentation states that the general assumptions follow Gadegaard et al. (2019).

- (c_{ij}): cost incurred when facility (i) serves customer (j).
- (f_i): fixed cost of opening facility (i).
- (d_j): demand of customer (j), drawn from (U(5,10)).

All parameters are non-negative integers.

The instance files also contain facility-capacity values. The model documentation explicitly says to ignore these values for this tri-objective formulation; no capacity constraint appears in the formulation reproduced above.

## Benchmark instances

The dataset contains **120 instances** in 12 subclasses, with 10 instances per subclass.

- Number of facilities: (5,10,ldots,60), in increments of five.
- Number of customers: (10,20,ldots,120), in increments of ten.
- Each instance has twice as many customers as facilities.

The filename convention is:

`<number_of_customers>_<number_of_facilities>_instance<number>`

For example, `10_5_instance1.txt` denotes an instance with 10 customers and 5 facilities.

## Instance file format

The original documentation specifies the following sequence:

1. Number of facilities and number of customers: (i), (j).
2. Facility fixed costs (f_i), followed by facility capacities (ignore the capacity values).
3. Customer demands (d_j).
4. Assignment costs (c_{ij}).

The files have the `.txt` extension. The original description does not specify the precise row/column traversal order of the assignment-cost matrix; the parser should be checked against the supplied files before assuming one.

## References

An, D., Parragh, S. N., Sinnl, M. and Tricoire, F. *Tri-objective facility location problem model*. Original accompanying model documentation.

Gadegaard, N. et al. (2019). Bi-objective single-source capacitated facility-location model. See the original model documentation's bibliography for the full bibliographic entry.

The TOFLP benchmark instances are included in the accompanying `TOFLP.zip` archive.
