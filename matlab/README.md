### Dependencies:
MATLAB R2019b or later is required. The `name=value` syntax in the examples requires R2021a or later; earlier supported releases can use `'name', value` pairs.

The MATLAB version of the repository makes use of the Multi-Parametric Toolbox 3.0 to handle projections of polytopes:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .

### Quick-start:
For more information on the quantities below please refer to [ALOT24](https://arxiv.org/abs/2107.08566). The main wrapper function takes nine positional matrix arguments followed by name-value options: `computeRCIS(A, B, E, Gx, Fx, Gu, Fu, Gw, Fw, name=value, ...)`.
  * `A`, `B`, and `E` are matrices defining the discrete-time linear system: `x^+ = Ax + Bu + Ew`.
  * `Gx` is a matrix and `Fx` is a vector defining the polyhedral safe set `Sx = {x \in \R^n | Gx x <= Fx}`.
  * `Gu` and `Fu` define input constraints `Su = {u \in \R^m | Gu u <= Fu}`. If no costraints use `Gu = []`, `Fu = []`.
  * `Gw` and `Fw` define the disturbance set `Sw = {w \in \R^k | Gw w <= Fw}`. If no disturbance use `E = []`, `Gw = []`, and `Fw = []`.
  * `lambda`: positive integer loop length; default `[]` (not specified).
  * `tau`: nonnegative integer transient length; default `0`.
  * `hierarchy_level`: positive integer hierarchy level; default `[]` (not specified).
  * `is_implicit`: logical scalar or numeric `0`/`1`; default `true`. Use `false` for explicit output via projection.

A nonempty `hierarchy_level` takes precedence over `lambda` and `tau`, whose values are then ignored without validation. Otherwise, a nonempty `lambda` is required. Empty `lambda` and `hierarchy_level` mean "not specified"; an explicitly empty `tau` is invalid when computing a single component.

The sequence length is `q = hierarchy_level` in hierarchy mode and `q = tau + lambda` in single-component mode. Hierarchy mode returns the individual components for `lambda = 1, ..., q`, with `tau = q - lambda`; it does not assemble their union.

```matlab
data = {A, B, E, Gx, Fx, Gu, Fu, Gw, Fw};

% Single component with tau = 0 and implicit output.
RCIS = computeRCIS(data{:}, lambda=3);

% Single component with a transient.
RCIS = computeRCIS(data{:}, lambda=3, tau=2);

% Hierarchy level 5 with explicit output.
RCIS = computeRCIS(data{:}, hierarchy_level=5, is_implicit=false);

% Hierarchy level 4 overrides lambda and tau.
RCIS = computeRCIS(data{:}, lambda=3, tau=2, hierarchy_level=4);

% Equivalent name-value syntax supported since R2019b.
RCIS = computeRCIS(data{:}, 'lambda', 3, 'tau', 2);
```

The previous positional configuration arguments are no longer supported; migrate calls to name-value options.

Output is a Polyhedron object, or an array of objects in hierarchy mode. Explicit output is in the original state space. Implicit output is lifted; its dimension depends on `q` and whether input constraints require state extension.

### Tests:
From the repository root in MATLAB, run:

```matlab
results = runtests('matlab/test');
assertSuccess(results);
```

Option-validation tests do not require MPT. Integration tests require MPT3 and Control System Toolbox on the MATLAB path; they are marked incomplete when those dependencies are absent. An incomplete run is not full validation of the API.
