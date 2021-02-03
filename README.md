# Controlled Invariant Sets in 2 Moves (cis2m)

This repository is dedicated to the computation of Controlled Invariant Sets (CISs) of Polyhedral Safe Sets for the class of Controllable Discrete-time Time-invariant Linear Systems. 

The proposed algorithm is based on an idea that works in two moves:
1. The problem is lifted in a higher dimensional space, where the Maximal Controlled Invariant Set (MCIS) is computed exactly and in closed-form.
2. The above MCIS is projected back to the original space, where it constitutes a CIS for the original problem.

Moreover, a hierarchy of CISs is established, which is parameterized by a positive integer denoting the level of the hierarchy. Typically higher levels of the hierarchy correspond to larger CISs (see Section 4.6 [AT20](https://dl.acm.org/doi/abs/10.1145/3365365.3382205) for more information).

### Related publications
1. Tzanis Anevlavis and Paulo Tabuada, 
"Computing controlled invariant sets in two moves", 
In 2019 IEEE Conference on Decision and Control (CDC). [AT19](https://ieeexplore.ieee.org/document/9029610). [Preprint](http://sites.google.com/a/g.ucla.edu/tzanis/home/anevlavisCDC2019.pdf).

2. Tzanis Anevlavis and Paulo Tabuada, 
"A simple hierarchy for computing controlled invariant sets", 
In Proceedings of the 23rd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'20). [AT20](https://dl.acm.org/doi/abs/10.1145/3365365.3382205).

3. T. Anevlavis, Z. Liu, N. Ozay and P. Tabuada,
"An enhanced hierarchy for (robust) controlled invariance",
In 2021 American Control Conference (ACC). (Accepted)

#####  More and cooler examples coming soon! Stay tuned!

For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.

### Dependencies:
The current version of the repository makes use of the Multi-Parametric Toolbox 3.0 to handle projections of polytopes:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .

### Quick-start:
The main wrapper function for this code is `computeCIS(A,B,Gx,Fx,L,Gu,Fu,E,Gw,Fw,method,verbose)`:
  * `A`, `B`, and `E` are matrices that define the discrete-time linear system `x^+ = Ax + Bu + Ew`.
  * `Gx` is a matrix, and `Fx` is a vector that define the polyhedral safe set `Sx = {x \in \R^n | Gx x <= Fx}`.
  * `L` is a positive integer denoting the level of the hierarchy. Typically larger values of `L` correspond to larger CISs (see Section 4.6 [AT20](https://dl.acm.org/doi/abs/10.1145/3365365.3382205) for more information).
  * `Gu` and `Fu` similarly define polyhedral input constraints `Su = {u \in \R^m | Gu u <= Fu}`. For unconstrained input use `Gu = []`, `Fu = []`.
  * `Gw` and `Fw` similarly define the disturbance set `Sw = {w \in \R^k | Gw w <= Fw}`. In absense of disturbance use `E = []`, `Gw = []`, `Fw = []`.
  * `method`: if not specified the default (best for most scenarios) algorithm will be selected (recommended). Legacy options for reference to specific paper results are: `CDC19`, `HSCC20`, `ACC21a`, `ACC21b`.
  * `verbose` is `0` for silent output, or `1`for verbose. Default is `0`.

The function will return a CIS of `Sx` in the form of a tuple `(cisA, cisb)` such that `CIS = {x \in \R^n | cisA x <= cisb}`.

### Citations:
If you used our algorithm for computing controlled invariant sets please cite as:
```latex
@inproceedings{anevlavis2020hscc
 author = {Anevlavis, Tzanis and Tabuada, Paulo},
 title = {A simple hierarchy for computing controlled invariant sets},
 booktitle = {Proceedings of the 23Rd ACM International Conference on Hybrid Systems: Computation and Control},
 series = {HSCC '20},
 year = {2020},
 isbn = {978-1-4503-7018-9/20/04},
 location = {Sydney, NSW, Australia},
 pages = {},
 numpages = {11},
 url = {http://doi.acm.org/10.1145/3365365.3382205},
 doi = {10.1145/3365365.3382205},
 acmid = {},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {Controller Synthesis, Safety, Controlled Invariance},
} 
```
or
```latex
@inproceedings{anevlavis2019cdc, 
author={T. {Anevlavis} and P. {Tabuada}}, 
booktitle={2019 IEEE 58th Conference on Decision and Control (CDC)}, 
title={Computing controlled invariant sets in two moves}, 
year={2019}, 
volume={}, 
number={}, 
doi = {10.1109/CDC40024.2019.9029610},
pages={6248-6254},}
```
