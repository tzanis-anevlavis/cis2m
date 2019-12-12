# cis2m
## Computing Controlled Invariant Sets in Two Moves

The provided code computes a Controlled Invariant Set for Discrete-time Linear Systems. It is an implementation in MATLAB of the algorithm proposed in:

Tzanis Anevlavis and Paulo Tabuada, 
"Computing controlled invariant sets in two moves", 
In 2019 IEEE Conference on Decision and Control (CDC). [Preprint](http://sites.google.com/a/g.ucla.edu/tzanis/home/anevlavisCDC2019.pdf).

In the `paper-examples` folder you can find the files that replicate the two examples of the above paper. 

#####  More and cooler examples coming soon! Stay tuned!

For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.

### Dependencies:
In this version of the code we make use of the Multi-Parametric Toolbox 3.0:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .

### Quick-start:
The main wrapper function for this code is `computeCIS(A,B,G,F,Gu,Fu,verbose)`. Simply provide as arguments the matrices `A` and `B` that define your discrete-time linear system `x^+ = Ax + Bu`, and `G` and `f` that define the polyhedral safe set `D = {x \in \R^n | Gx <= F}`. The function will return a controlled invariant subset of D in the form of a matrix `[Gcis Fcis]` such that `CIS = {x \in \R^n | Gcis x <= Fcis}`. Matrices `Gu` and `Fu` define similarly polyhedral constraints on the input `u`, i.e., `Du = {u \in \R | Gu u <= Fu}`. For unconstrained input use `[]` at each entry `Gu`, `Fu`.

A more detailed read-me will follow..
