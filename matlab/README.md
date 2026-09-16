## API
For MATLAB dependencies see the Dependencies section at the bottom.
### Quick-start
The main wrapper function takes seven positional matrix arguments followed by name-value options: `computeRCIS(A, B, E, Gxu, Fxu, Gw, Fw, name=value, ...)`:
  * `A`, `B`, and `E` are matrices defining the discrete-time linear system: `x^+ = Ax + Bu + Ew`.
  * `Gxu` and the column vector `Fxu` define the joint safe set `Sxu = {(x,u) | Gxu*[x;u] <= Fxu}`. The first `n` columns multiply the state; the last `m` multiply the physical input. Rows may couple states and inputs.
  * `Gw` and `Fw` define the disturbance set `Sw = {w \in \R^k | Gw w <= Fw}`. If no disturbance use `E = []`, `Gw = []`, and `Fw = []`.
  * `lambda`: positive integer loop length; default `[]` (not specified).
  * `tau`: nonnegative integer transient length; default `0`.
  * `hierarchy_level`: positive integer hierarchy level; default `[]` (not specified).
  * `is_implicit`: logical scalar or numeric `0`/`1`; default `true`. Use `false` for explicit output via projection.

A nonempty `hierarchy_level` takes precedence over `lambda` and `tau`, whose values are then ignored without validation. Otherwise, a nonempty `lambda` is required. Empty `lambda` and `hierarchy_level` mean "not specified"; an explicitly empty `tau` is invalid when computing a single component.

The sequence length is `q = hierarchy_level` in hierarchy mode and `q = tau + lambda` in single-component mode. Hierarchy mode returns the individual components for `lambda = 1, ..., q`, with `tau = q - lambda`; it does not assemble their union.

For more information on the quantities above refer to [ALOT24](https://ieeexplore.ieee.org/document/10328804) / [arxiv](https://arxiv.org/abs/2107.08566).

Sample API usage:
```matlab
data = {A, B, E, Gxu, Fxu, Gw, Fw};
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

The API supports joint state-input constraints. To use separate constraints:

```matlab
Gxu = blkdiag(Gx, Gu);
Fxu = [Fx; Fu];
% State-only constraints, with m = size(B, 2).
Gxu = [Gx zeros(size(Gx, 1), size(B, 2))];
Fxu = Fx;
```

Output is a Polyhedron object, or an array of objects in hierarchy mode. Explicit output is in the original state space. Implicit output uses `[x; v]` coordinates and has dimension `n + m * q`. The virtual inputs `v` are grouped by input channel, with `q` entries per channel, and `r = H*v`.

The optional second output `A_lifted` contains the nominal transition dynamics in `[x; v]` coordinates for either output mode. For a single `(tau, lambda)` component it is one sparse matrix. In hierarchy mode it is a cell array: `A_lifted{i}` corresponds to `RCIS(i)`, with `lambda = i` and `tau = q - i`.

### Lifted coordinates and the physical input

The implicit RCIS is a set in `[x; v]`, i.e., the (state, virtual input)-space, and not in the  `[x; u]`, i.e., the (state, physical input)-space.

The vector `v` is an auxiliary controller state, part of the lifted state and system:
$$
\begin{bmatrix}z^+\\v^+\end{bmatrix}
=
\begin{bmatrix}A&BH\\0&P\end{bmatrix}
\begin{bmatrix}z\\v\end{bmatrix}.
$$
The physical input `u` is an output of the lifted state and system:
$$
u(x,v)=B_m^{-1}(Hv-A_mTx).
$$
Note, $v\in\mathbb R^{mq}$, but $u\in\mathbb R^m$, so mapping $v$ to $u$ discards the future lasso samples and is not an invertible coordinate transformation.

The transformed lifted dynamics are:
$$
\begin{bmatrix}x^+\\v^+\end{bmatrix}
=
A_{xv}
\begin{bmatrix}x\\v\end{bmatrix}
=
\begin{bmatrix}
A-BB_m^{-1}A_mT & BB_m^{-1}H\\
0&P
\end{bmatrix}
\begin{bmatrix}x\\v\end{bmatrix}
$$
which gives:
$$
x^+=Ax+Bu(x,v).
$$

If an $(x,u)$ is desired, it would be the linear image:
$$
\begin{bmatrix}x\\u\end{bmatrix}
=
\begin{bmatrix}
I&0\\
-B_m^{-1}A_mT&B_m^{-1}H
\end{bmatrix}
\begin{bmatrix}x\\v\end{bmatrix}.
$$
For simplicitly in the above, the external disturbance $w$ and matrix $E$ are omitted.

For a feasible lifted point `[x; v]`, recover the physical input using:
```matlab
[~, ~, ~, ~, ~, T, ~, Am, Bm] = ...
    transformToBrunovskyNormalForm(A, B, E, Gxu, Fxu);
H = kron(speye(size(B, 2)), [1 sparse(1, q - 1)]);
u = Bm \ (H*v - Am*T*x);
```

----
#### Supervisory control example
The implicit RCIS can be used for supervisory control in two related ways.
##### 1. One-step supervision using the implicit set
Select a safe physical input and certify that the successor state belongs to
$\operatorname{proj}_x(\mathcal C_{xv})$:
$$
\begin{aligned}
\min_{u,v^+}\quad &\|u-u_{\mathrm{ff}}\|\\
\text{s.t.}\quad &[x;u]\in S_{xu},\\
&[Ax+Bu;v^+]\in\mathcal C_{xv}.
\end{aligned}
$$
Here $v^+$ is only a witness that the successor state belongs to the
projection of the implicit RCIS. The current input $u$ is selected
independently of the lifted controller.

##### 2. Supervision using the lifted controller directly
Alternatively, select a current virtual controller state $v$ and obtain the
physical input through the feedback transformation:
$$
\begin{aligned}
\min_{u,v}\quad &\|u-u_{\mathrm{ff}}\|\\
\text{s.t.}\quad &[x;u]\in S_{xu},\\
&[x;v]\in\mathcal C_{xv},\\
&B_m u + A_mTx = Hv.
\end{aligned}
$$
The constraint $[x;u]\in S_{xu}$ is already implied by the last two
constraints and the construction of $\mathcal C_{xv}$, but is shown explicitly
for clarity. Invariance then guarantees
$$
[Ax+Bu;Pv]\in\mathcal C_{xv}.
$$
Eliminating $u$ gives the equivalent formulation:
$$
\begin{aligned}
\min_v\quad&
\left\|B_m^{-1}(Hv-A_mTx)-u_{\mathrm{ff}}\right\|\\
\text{s.t.}\quad &[x;v]\in\mathcal C_{xv}.
\end{aligned}
$$

## Tests
From the repository root in MATLAB, run:
```matlab
results = runtests('matlab/test');
assertSuccess(results);
```
Option validation and lifted-matrix tests do not require MPT. Transformation tests require Control System Toolbox. Disturbance and integration tests require MPT3; end-to-end tests also need Control System Toolbox. Tests are marked incomplete when their dependencies are absent. An incomplete run is not full validation of the API.

## Dependencies
MATLAB R2019b or later is required. The `name=value` syntax in the examples requires R2021a or later; earlier supported releases can use `'name', value` pairs.

The MATLAB version of the repository makes use of the Multi-Parametric Toolbox 3.0 to handle projections of polytopes:
M. Herceg, M. Kvasnica, C. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510. http://control.ee.ethz.ch/mpt .