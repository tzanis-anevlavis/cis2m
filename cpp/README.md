# C++ implementation

This directory contains the C++ implementation of `cis2m`.
See the [repository README](../README.md) for the paper and project
overview.

## Getting started

### Dependencies

- A C++14 compiler and CMake 3.15 or newer
- Eigen3 and OR-Tools for the library
- GoogleTest 1.14 or newer for tests. If it is not installed, CMake downloads
  GoogleTest 1.14 into the build tree; it is not installed with cis2m.

The build script below uses Homebrew, so it requires Homebrew and, on macOS,
Apple Command Line Tools. A manual CMake build can use dependencies installed
through another package manager.

### Build

From the repository root, the Homebrew-based build is:

```sh
cd cpp
./build_cis2m.sh
```

The script updates Homebrew, installs missing build dependencies, locates a
compatible OR-Tools installation or builds one, then configures and builds
cis2m. **It recreates `cpp/build` on every run.** It does not run the tests or
install cis2m.

To use dependencies that are already installed and discoverable by CMake:

```sh
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build --parallel
```

Pass `-DCMAKE_PREFIX_PATH=/path/to/dependencies` at configuration time if
Eigen3 or OR-Tools is installed outside CMake's default search paths. To build
only the library, also pass `-DBUILD_TESTING=OFF`.

### Tests

After building with testing enabled, run from `cpp/build`:

```sh
ctest --output-on-failure
```

The suite includes polyhedron, Brunovsky transformation, CIS generator, and
installed-package tests. To run the polyhedron and transformation tests with
the installed-package check separately:

```sh
ctest --output-on-failure -R 'test_(hpolyhedron|brunovskytransformation|cmake_package)'
```

The package test installs into a temporary directory under the build tree,
then configures and runs an independent CMake consumer.

The CIS generator tests include dense, nontrivial joint state-input safe sets
under dense coordinate and feedback transformations. Nominal CIS and robust
RCIS invariance are verified for single-input systems with up to 20 states and
a controllability index of 20, and for two-input systems with up to 20 states.
These tests exercise high-dimensional implicit-set construction without an
explicit projection step.

### Install

To install the built library, public headers, and CMake package files to a
chosen prefix:

```sh
cmake --install cpp/build --prefix "$HOME/local/cis2m"
```

A downstream CMake project can use `find_package(cis2m REQUIRED CONFIG)` and
link `cis2m::cis2m`. Add the installation prefix to `CMAKE_PREFIX_PATH` if
needed. Eigen3 is a public dependency; the OR-Tools library must also be
available at runtime.

## Controlled invariant set generator API

`ControlledInvariantSetGenerator` computes invariant sets for:
$$
x^+ = A x + B u + E w.
$$

Construct it with `(A, B)` for a nominal system or `(A, B, E)` when the model
accepts disturbances. The pair `(A, B)` must be controllable, and `B` must
have full column rank.

The safe set can be supplied as an `HPolyhedron` over either `x` or `[x; u]`.
State-only constraints leave `u` unconstrained. A nonempty disturbance set
must be bounded and have dimension `E.cols()`. Matrix overloads also accept
the inequalities `Gxu * [x; u] <= Fxu` and `Gw * w <= Fw` directly.

### Options and results

Set `CISOptions::tau` and `CISOptions::lambda` to compute one lasso component,
where `lambda > 0`. Alternatively, set `hierarchy_level = q` to compute all
`q` components with `lambda = 1, ..., q` and `tau = q - lambda`; explicitly
provided `tau` and `lambda` are then ignored.

The default `is_implicit = true` returns the set in lifted `[x; v]`
coordinates, with `v` in $\mathbb{R}^{m(\tau+\lambda)}$. Each `CISComponent`
also provides

```text
[x+; v+] = lifted_dynamics * [x; v] + lifted_disturbance * w
         u = input_from_state * x + input_from_virtual * v.
```

Setting `is_implicit = false` projects the set onto `x`. Projection uses
Fourier-Motzkin elimination and can be substantially more expensive and
numerically sensitive than returning the implicit set.

### Example

```cpp
#include <cis2m/cis_generator.hpp>

#include <Eigen/Dense>

using cis2m::ControlledInvariantSetGenerator;
using cis2m::CISOptions;
using cis2m::HPolyhedron;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
    MatrixXd A(2, 2);
    A << 1.0, 1.0,
         0.0, 1.0;
    MatrixXd B(2, 1);
    B << 0.0, 1.0;
    MatrixXd E(2, 1);
    E << 0.0, 0.1;

    // Joint box constraints on [x; u].
    MatrixXd Gxu(6, 3);
    Gxu <<  1.0,  0.0,  0.0,
           -1.0,  0.0,  0.0,
            0.0,  1.0,  0.0,
            0.0, -1.0,  0.0,
            0.0,  0.0,  1.0,
            0.0,  0.0, -1.0;
    VectorXd Fxu = VectorXd::Ones(6);
    HPolyhedron safe_set(Gxu, Fxu);

    // Disturbance interval |w| <= 0.05.
    MatrixXd Gw(2, 1);
    Gw << 1.0, -1.0;
    VectorXd Fw(2);
    Fw << 0.05, 0.05;
    HPolyhedron disturbance_set(Gw, Fw);

    CISOptions options;
    options.tau = 1;
    options.lambda = 2;
    options.is_implicit = true;

    ControlledInvariantSetGenerator generator(A, B, E);
    const auto components = generator.Compute(safe_set, disturbance_set, options);
    const auto& component = components.front();

    // component.set constrains [x; v]. The corresponding physical input is:
    VectorXd x = VectorXd::Zero(2);
    VectorXd v = VectorXd::Zero(3);
    VectorXd u = component.input_from_state * x +
                 component.input_from_virtual * v;
    return u.allFinite() ? 0 : 1;
}
```

Omit `E` and the disturbance set for nominal computation:

```cpp
ControlledInvariantSetGenerator nominal_generator(A, B);
const auto nominal_components = nominal_generator.Compute(safe_set, options);
```

Invalid systems, incompatible sets, invalid options, and detected numerical
failures are reported with exceptions.
