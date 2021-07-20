# Controlled Invariant Sets in 2 Moves (cis2m) and closed-form implicit representations

This repository is dedicated to the computation of (Robust) Controlled Invariant Sets (RCISs) of Polyhedral Safe Sets for the class of Controllable Discrete-time Time-invariant Linear Systems. 

The proposed algorithm operates in two modes:
1. it returns a *closed-form* expression for an implicit (higher dimensional) representation of an RCIS in the space of states and finite input sequences. 
2. it optionally returns an explicit RCIS by projecting the implicit representation to the original state space.

The closed-form expression allows for scalability and is suitable for online applications. Moreover, a hierarchy of RCISs is established, which is parameterized by a positive integer (level of the hierarchy).

The following video demonstrates how our approach guarantees collision-free trajectories when supervising a Crazyflie 2.0 quadrotor for the task of obstacle avoidance!

![drone_supervision](https://user-images.githubusercontent.com/26322321/110282721-d4150700-7f93-11eb-8537-2edab340b7ac.gif)

Full video [here](https://tinyurl.com/drone-supervision-cis).


### Related publications
1. T. Anevlavis, Z. Liu, N. Ozay and P. Tabuada,
"Controlled invariant sets: implicit closed-form representations and applications",
arXiv:2107.08566 [math.OC], 2021. [ALOT21](https://arxiv.org/abs/2107.08566)

2. Tzanis Anevlavis and Paulo Tabuada, 
"Computing controlled invariant sets in two moves", 
In 2019 IEEE Conference on Decision and Control (CDC). [AT19](https://ieeexplore.ieee.org/document/9029610).

3. Tzanis Anevlavis and Paulo Tabuada, 
"A simple hierarchy for computing controlled invariant sets", 
In Proceedings of the 23rd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'20). [AT20](https://dl.acm.org/doi/abs/10.1145/3365365.3382205).

4. T. Anevlavis, Z. Liu, N. Ozay and P. Tabuada,
"An enhanced hierarchy for (robust) controlled invariance",
In 2021 American Control Conference (ACC). (Accepted)

5. Luigi Pannocchi, Tzanis Anevlavis, and Paulo Tabuada,
"Trust your supervisor: quadrotor obstacle avoidance using controlled invariant sets", 
In 2021 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2021. (Accepted)

#####  More and cooler examples coming soon! Stay tuned!

For any comments contact Tzanis Anevlavis @ t.anevlavis@ucla.edu.

### Dependencies:
The MATLAB version of the repository makes use of the Multi-Parametric Toolbox 3.0 to handle projections of polytopes:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .

### Quick-start:
For more information on the quantities below please refer to  [ALOT21](https://arxiv.org/abs/2107.08566). The main wrapper function for this code is `computeRCIS(A,B,Gx,Fx,Gu,Fu,E,Gw,Fw,implicit,L,T)`:
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

### Citations:
If you used our algorithm for computing controlled invariant sets please cite as:
```latex
@misc{anevlavis2021controlled,
      title={Controlled invariant sets: implicit closed-form representations and applications}, 
      author={Tzanis Anevlavis and Zexiang Liu and Necmiye Ozay and Paulo Tabuada},
      year={2021},
      eprint={2107.08566},
      archivePrefix={arXiv},
      primaryClass={math.OC}
}
```