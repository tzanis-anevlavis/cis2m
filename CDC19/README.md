## Computing Controlled Invariant Sets in Two Moves

The provided code computes a Controlled Invariant Set (CIS) of a Polyhedron for the class of Controllable Discrete-time Linear Systems. It is an implementation in MATLAB of the algorithm proposed in:

Tzanis Anevlavis and Paulo Tabuada, 
"Computing controlled invariant sets in two moves", 
In 2019 IEEE Conference on Decision and Control (CDC). [Preprint](http://sites.google.com/a/g.ucla.edu/tzanis/home/anevlavisCDC2019.pdf).

In the `paper-examples` folder you can find the files that replicate the two examples of the above paper. 

#####  More and cooler examples coming soon! Stay tuned!

For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.

### Dependencies:
In this version of the code we make use of the Multi-Parametric Toolbox 3.0:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .

### Related publications
1. Tzanis Anevlavis and Paulo Tabuada, 
"Computing controlled invariant sets in two moves", 
In 2019 IEEE Conference on Decision and Control (CDC). [Preprint](http://sites.google.com/a/g.ucla.edu/tzanis/home/anevlavisCDC2019.pdf).

2. Tzanis Anevlavis and Paulo Tabuada, 
"A simple hierarchy for computing controlled invariant sets", 
In Proceedings of the 23rd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'20).

#####  More and cooler examples coming soon! Stay tuned!

For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.

### Dependencies:
The current version of the repository makes use of the Multi-Parametric Toolbox 3.0:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .

### Quick-start:
The main wrapper function for this code is `computeCIS(A,B,G,F,Gu,Fu,method,verbose)`:
  * `A` and `B` are matrices that define your discrete-time linear system `x^+ = Ax + Bu`.
  * `G` and `f` are matrices that define the polyhedral safe set `D = {x \in \R^n | Gx <= F}`.
  * `Gu` and `fu` are similarly matrices that define the polyhedral constraints on the input `u`, i.e., `Du = {u \in \R | Gu u <= Fu}`. For unconstrained input use `[]` at each entry `Gu`, `Fu`.
  * `method` is either `'CDC19'`, or `'HSCC20'`. More details for each method in the corresponding folders.
  * `verbose` is `0` for silent output, or `1`for verbose. Default is `0`.

The function will return a controlled invariant subset of D in the form of a matrix `[Gcis Fcis]` such that `CIS = {x \in \R^n | Gcis x <= Fcis}`.

### Citations:
If you used this algorithm for computing controlled invariant sets please cite as:
```latex
@inproceedings{AT2019cis2m,
author={T. Anevlavis and P. Tabuada}, 
booktitle={2019 IEEE Conference on Decision and Control (CDC)}, 
title={Computing controlled invariance sets in two moves}, 
year={2019}, 
volume={}, 
number={}, 
pages={6249-6254}, 
keywords={}, 
doi={}, 
url={},
ISSN={}, 
month={Dec},} 
```

