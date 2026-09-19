#!/usr/bin/env bash

set -e  # Exit immediately on error
set -o pipefail

echo "==> Ensuring build dependencies are installed via Homebrew..."
brew update

brew install cmake git pkgconf || true

echo "==> Installing Eigen3..."
brew install eigen || true

# === CONFIGURATION ===
ORTOOLS_VERSION="9.14"

# Get script directory and derive paths from it
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CIS2M_ROOT="$(dirname "$SCRIPT_DIR")"
PARENT_DIR="$(dirname "$CIS2M_ROOT")"
CIS2M_SRC_DIR="$SCRIPT_DIR"
CIS2M_BUILD_DIR="$CIS2M_SRC_DIR/build"
HOST_ARCH="$(uname -m)"
HOST_OS="$(uname -s)"
MACOS_SDK_PATH=""
CMAKE_DARWIN_FLAGS=()

if [[ "$HOST_OS" == "Darwin" ]] && command -v xcrun >/dev/null 2>&1; then
  MACOS_SDK_PATH="$(xcrun --sdk macosx --show-sdk-path 2>/dev/null || true)"
  if [[ -n "$MACOS_SDK_PATH" ]]; then
    CMAKE_DARWIN_FLAGS=(-DCMAKE_OSX_SYSROOT="$MACOS_SDK_PATH")
  fi
fi

BREW_PREFIX="$(brew --prefix)"
BREW_ORTOOLS_PREFIX="$(brew --prefix or-tools 2>/dev/null || true)"
ORTOOLS_SRC_DIR="$PARENT_DIR/ortools"
ORTOOLS_FALLBACK_INSTALL_DIR="$PARENT_DIR/local/ortools"

has_ortools_layout() {
  local candidate="$1"
  [[ -f "$candidate/lib/cmake/ortools/ortoolsConfig.cmake" ]] && \
  [[ -f "$candidate/include/ortools/linear_solver/linear_solver.h" ]]
}

ortools_library_for_arch_check() {
  local candidate="$1"
  local lib
  for lib in \
    "$candidate/lib/libortools.dylib" \
    "$candidate/lib/libortools.9.dylib" \
    "$candidate/lib/libortools.9.0.dylib" \
    "$candidate/lib/libortools.a"; do
    if [[ -f "$lib" ]]; then
      echo "$lib"
      return 0
    fi
  done
  return 1
}

is_ortools_arch_compatible() {
  local candidate="$1"

  # Non-macOS: skip architecture filtering here.
  if [[ "$HOST_OS" != "Darwin" ]]; then
    return 0
  fi

  local lib_path
  if ! lib_path="$(ortools_library_for_arch_check "$candidate")"; then
    echo "==> Warning: Could not find OR-Tools core library under $candidate/lib for architecture check." >&2
    return 0
  fi

  local lib_arches
  lib_arches="$(lipo -archs "$lib_path" 2>/dev/null || true)"
  if [[ -z "$lib_arches" ]]; then
    echo "==> Warning: Could not determine library architectures for $lib_path." >&2
    return 0
  fi

  if [[ " $lib_arches " == *" $HOST_ARCH "* ]]; then
    return 0
  fi

  echo "==> Skipping OR-Tools at $candidate (library arch: $lib_arches, host arch: $HOST_ARCH)" >&2
  return 1
}

validate_ortools_source_checkout() {
  local expected_tag="v$ORTOOLS_VERSION"
  local actual_tag

  if ! git -C "$ORTOOLS_SRC_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "Error: Existing OR-Tools source directory is not a Git checkout: $ORTOOLS_SRC_DIR" >&2
    echo "Remove or relocate it so the build script can clone $expected_tag." >&2
    return 1
  fi

  actual_tag="$(git -C "$ORTOOLS_SRC_DIR" describe --tags --exact-match HEAD 2>/dev/null || true)"
  if [[ "$actual_tag" != "$expected_tag" ]]; then
    if [[ -z "$actual_tag" ]]; then
      actual_tag="untagged commit $(git -C "$ORTOOLS_SRC_DIR" rev-parse --short HEAD)"
    fi
    echo "Error: OR-Tools checkout at $ORTOOLS_SRC_DIR is $actual_tag; expected $expected_tag." >&2
    echo "Remove or relocate it so the build script can clone the requested version." >&2
    return 1
  fi
}

find_compatible_ortools_install() {
  local candidate
  local -a candidates=(
    "$HOME/local/ortools"
    "$PARENT_DIR/local/ortools"
    "$BREW_ORTOOLS_PREFIX"
    "$BREW_PREFIX/opt/or-tools"
    "/opt/homebrew/opt/or-tools"
    "/opt/homebrew"
    "/usr/local"
  )

  for candidate in "${candidates[@]}"; do
    [[ -z "$candidate" ]] && continue
    if has_ortools_layout "$candidate"; then
      if is_ortools_arch_compatible "$candidate"; then
        echo "$candidate"
        return 0
      fi
    fi
  done

  return 1
}

# Try to locate existing OR-Tools installation
ORTOOLS_INSTALL_DIR=""
if ORTOOLS_INSTALL_DIR="$(find_compatible_ortools_install)"; then
  echo "==> Found compatible OR-Tools at $ORTOOLS_INSTALL_DIR"
else
  echo "==> No compatible OR-Tools installation found. Trying Homebrew package..."
  brew install or-tools || true
  BREW_ORTOOLS_PREFIX="$(brew --prefix or-tools 2>/dev/null || true)"

  if ORTOOLS_INSTALL_DIR="$(find_compatible_ortools_install)"; then
    echo "==> Found compatible OR-Tools at $ORTOOLS_INSTALL_DIR"
  else
    ORTOOLS_INSTALL_DIR="$ORTOOLS_FALLBACK_INSTALL_DIR"
    echo "==> Will build OR-Tools from source at $ORTOOLS_SRC_DIR"
  fi
fi

# Build OR-Tools if not already found
if has_ortools_layout "$ORTOOLS_INSTALL_DIR" && is_ortools_arch_compatible "$ORTOOLS_INSTALL_DIR"; then
  echo "==> Using OR-Tools installation at $ORTOOLS_INSTALL_DIR"
else
  echo "==> Building OR-Tools from source..."
  
  # Clone the requested release, or verify that an existing checkout matches it.
  if [[ ! -d "$ORTOOLS_SRC_DIR" ]]; then
    mkdir -p "$(dirname "$ORTOOLS_SRC_DIR")"
    git clone --branch "v$ORTOOLS_VERSION" --depth 1 \
      https://github.com/google/or-tools.git "$ORTOOLS_SRC_DIR"
  else
    validate_ortools_source_checkout
  fi
  
  cd "$ORTOOLS_SRC_DIR" || { echo "Error: Cannot change to OR-Tools directory $ORTOOLS_SRC_DIR"; exit 1; }
  rm -rf build
  cmake -S . -B build \
    -DBUILD_DEPS=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$ORTOOLS_INSTALL_DIR" \
    "${CMAKE_DARWIN_FLAGS[@]}"
  
  cmake --build build --parallel
  cmake --install build --prefix "$ORTOOLS_INSTALL_DIR"
  
  echo "==> OR-Tools built and installed to $ORTOOLS_INSTALL_DIR"
fi

# === Build cis2m ===
echo "==> Building cis2m..."

cd "$CIS2M_SRC_DIR" || { echo "Error: Cannot change to cis2m directory $CIS2M_SRC_DIR"; exit 1; }
rm -rf "$CIS2M_BUILD_DIR"

# Prefer Homebrew pkg-config metadata for COIN-OR dependencies when using
# Homebrew OR-Tools, and avoid mixing with legacy /usr/local CMake packages.
export PKG_CONFIG_PATH="$BREW_PREFIX/lib/pkgconfig:$BREW_PREFIX/share/pkgconfig:${PKG_CONFIG_PATH:-}"
if command -v pkg-config >/dev/null 2>&1; then
  PKG_CONFIG_BIN="$(command -v pkg-config)"
elif command -v pkgconf >/dev/null 2>&1; then
  PKG_CONFIG_BIN="$(command -v pkgconf)"
else
  PKG_CONFIG_BIN="$BREW_PREFIX/bin/pkgconf"
fi

cmake -S . -B "$CIS2M_BUILD_DIR" \
  -DCMAKE_PREFIX_PATH="$ORTOOLS_INSTALL_DIR;$BREW_PREFIX" \
  -Dortools_DIR="$ORTOOLS_INSTALL_DIR/lib/cmake/ortools" \
  -DPKG_CONFIG_EXECUTABLE="$PKG_CONFIG_BIN" \
  -DCbc_NO_Cbc_CMAKE=ON \
  -DClp_NO_Clp_CMAKE=ON \
  -DCMAKE_BUILD_RPATH="$ORTOOLS_INSTALL_DIR/lib" \
  -DCMAKE_INSTALL_RPATH="$ORTOOLS_INSTALL_DIR/lib" \
  "${CMAKE_DARWIN_FLAGS[@]}"

cmake --build "$CIS2M_BUILD_DIR" --parallel

echo "✅ cis2m built successfully."
