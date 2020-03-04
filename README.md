# Computing Controlled Invariant Sets in Two Moves (cis2m)

This repository contains different algorithms that compute Controlled Invariant Sets (CISs) of Polyhedra for the class of Controllable Discrete-time Time-invariant Linear Systems. 

The proposed methods are based on an idea that works in two moves:
1. The problem is lifted in a higher dimensional space, where the Maximal Controlled Invariant Set (MCIS) is computed exactly and in closed-form.
2. The above MCIS is projected back to the original space, where it constitutes a CIS for the original problem.

The following different methods are available to choose from:
1. [CDC19](https://github.com/janis10/cis2m/tree/master/paper-archive/CDC19) - typically efficient for systems up to 7-8 dimensions.
2. [HSCC20](https://github.com/janis10/cis2m/tree/master/paper-archive/HSCC20) - a hierarchy for computing larger controlled invariant sets. 
3. [CDC20a] - same hierarchy as HSCC20, with superior computational efficiency. 
4. [CDC20b] - a less conservative hierarchy guaranteed to always include the results of [HSCC20/CDC20a] (preferred method).

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
  * `method` is either `'CDC19'`, `'HSCC20'`, `'CDC20a'`, or `'CDC20b'`. More details for each method in the corresponding folders.
  * `verbose` is `0` for silent output, or `1`for verbose. Default is `0`.

The function will return a controlled invariant subset of D in the form of a matrix `[Gcis Fcis]` such that `CIS = {x \in \R^n | Gcis x <= Fcis}`.

### Citations:
If you used our algorithm for computing controlled invariant sets please cite as:
```latex
@inproceedings{AT2020cis2mHier
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
depending on the version you used!
