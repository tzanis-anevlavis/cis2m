# Controlled invariant sets in 2 moves (`cis2m`)

`cis2m` computes controlled invariant sets (CISs) and robust controlled
invariant sets (RCISs) for controllable, discrete-time, time-invariant linear
systems with polyhedral state or joint state-input safety constraints.

The method considers systems of the form
$$
x^+ = A x + B u + E w,
$$
and provides two outputs:

1. A closed-form implicit invariant set $\mathcal C_{xv}$ in the original state
   $x$ and a finite-dimensional virtual input $v$. The variable $v$ encodes an
   eventually periodic virtual-input sequence.
2. An optional explicit CIS or RCIS obtained by projecting $\mathcal C_{xv}$
   onto the original state space.

The implicit construction avoids fixed-point set iterations and supports a
hierarchy of invariant sets indexed by a positive integer. Increasing the
hierarchy level can produce larger projected sets at the cost of a larger
lifted space.

## Demonstration

The following experiment uses an implicit RCIS to supervise a Crazyflie 2.0
quadrotor and guarantee collision-free obstacle avoidance.

![Crazyflie obstacle-avoidance supervision](https://user-images.githubusercontent.com/26322321/110282721-d4150700-7f93-11eb-8537-2edab340b7ac.gif)

Watch the [full demonstration](https://tinyurl.com/drone-supervision-cis).

## Implementations

- [MATLAB](./matlab/README.md) provides the reference implementation, API
  documentation, and MATLAB/MPT3 test instructions.
- [C++](./cpp/README.md) provides the C++ library, build and installation
  instructions, API documentation, and examples.

The C++ test suite verifies positive invariance of nominal implicit CISs and
robust positive invariance of implicit RCISs for dense, nontrivial joint
state-input constraints and dense coordinate and feedback transformations. Its
high-dimensional cases include single-input and two-input systems with up to
20 states, including a single-input system with controllability index 20.

Explicit projection is optional. It uses Fourier-Motzkin elimination and can
be substantially more expensive and numerically sensitive than retaining the
implicit representation, particularly in high dimensions.

## Getting Started

- For MATLAB, follow the [quick-start and dependency instructions](./matlab/README.md).
- For C++, follow the [build, test, and installation instructions](./cpp/README.md).

## Repository Structure

- [`cpp`](./cpp): C++ implementation and tests.
- [`matlab`](./matlab): MATLAB implementation and tests.
- [`paper-archive`](./paper-archive): files used to reproduce results from the
  associated publications. These files may not track later MATLAB and C++
  changes.

## Related Publications

1. T. Anevlavis, Z. Liu, N. Ozay, and P. Tabuada, "Controlled Invariant Sets:
   Implicit Closed-Form Representations and Applications," *IEEE Transactions
   on Automatic Control*, vol. 69, no. 7, pp. 4506-4521, 2024.
   [ALOT24](https://doi.org/10.1109/TAC.2023.3336819)
2. T. Anevlavis and P. Tabuada, "Computing Controlled Invariant Sets in Two
   Moves," *2019 IEEE Conference on Decision and Control (CDC)*, 2019.
   [AT19](https://ieeexplore.ieee.org/document/9029610)
3. T. Anevlavis and P. Tabuada, "A Simple Hierarchy for Computing Controlled
   Invariant Sets," *Proceedings of the 23rd ACM International Conference on
   Hybrid Systems: Computation and Control (HSCC)*, 2020.
   [AT20](https://doi.org/10.1145/3365365.3382205)
4. T. Anevlavis, Z. Liu, N. Ozay, and P. Tabuada, "An Enhanced Hierarchy for
   (Robust) Controlled Invariance," *2021 American Control Conference (ACC)*,
   pp. 4860-4865, 2021.
   [ALOT21](https://ieeexplore.ieee.org/document/9483217)
5. L. Pannocchi, T. Anevlavis, and P. Tabuada, "Trust Your Supervisor:
   Quadrotor Obstacle Avoidance Using Controlled Invariant Sets," *2021
   IEEE/RSJ International Conference on Intelligent Robots and Systems
   (IROS)*, pp. 9219-9224, 2021.
   [PAT21](https://doi.org/10.1109/IROS51168.2021.9636485)

## Citation

If you use this algorithm, please cite:

```bibtex
@article{anevlavis2024controlled,
  author  = {Tzanis Anevlavis and Zexiang Liu and Necmiye Ozay and Paulo Tabuada},
  title   = {Controlled Invariant Sets: Implicit Closed-Form Representations and Applications},
  journal = {IEEE Transactions on Automatic Control},
  year    = {2024},
  volume  = {69},
  number  = {7},
  pages   = {4506--4521},
  doi     = {10.1109/TAC.2023.3336819}
}
```

## Contact

For comments or questions, contact Tzanis Anevlavis at
`t.anevlavis@ucla.edu`.

## License

This project is licensed under the [GNU General Public License v3.0](./LICENSE).
