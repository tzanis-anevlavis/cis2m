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
installed-package tests. The CIS generator is under refactoring, so its test
result does not yet establish correctness of the complete method. To run the
finished components and the installed-package check separately:

```sh
ctest --output-on-failure -R 'test_(hpolyhedron|brunovskytransformation|cmake_package)'
```

The package test installs into a temporary directory under the build tree,
then configures and runs an independent CMake consumer.

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
