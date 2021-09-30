### Dependencies:
The MATLAB version of the repository makes use of the Multi-Parametric Toolbox 3.0 to handle projections of polytopes:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .

### Quick-start:
For more information on the quantities below please refer to  [ALOT21b](https://arxiv.org/abs/2107.08566). The main wrapper function for this code is `computeRCIS(A,B,Gx,Fx,Gu,Fu,E,Gw,Fw,implicit,L,T)`:
  * `A`, `B`, and `E` are matrices defining the discrete-time linear system: `x^+ = Ax + Bu + Ew`.
  * `Gx` is a matrix and `Fx` is a vector defining the polyhedral safe set `Sx = {x \in \R^n | Gx x <= Fx}`.
  * `Gu` and `Fu` define input constraints `Su = {u \in \R^m | Gu u <= Fu}`. If no costraints use `Gu = []`, `Fu = []`.
  * `Gw` and `Fw` define the disturbance set `Sw = {w \in \R^k | Gw w <= Fw}`. If no disturbance use `E = []`, `Gw = []`, and `Fw = []`.
  * `implicit \in {0,1}`.
    * If `implicit=0`, then `explicit RCIS = {x \in \R^n | rcisA x <= rcisb}`.
    * If `implicit=1`, then `implicit RCIS = {(x,u,v) \in \R^n x \R^m x \R^{m(T+L)} | rcisA(x,u,v) <= rcisb}`.
  * If only `L` is specified, then `L : L-th level of hierarchy`.
  * If both `L` and `T` are specified:
    * `L (lambda)`: loop of eventually periodic input sequence.
    * `T (tau)`:    transient of eventually periodic input sequence.

Output is a Polyhedron object. If `implicit=0`, then the output is an `explicit RCIS` as above. Else, it is an `implicit RCIS`.