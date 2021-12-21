### Dependencies:
The C++ version of the library makes use of [OR-Tools](https://developers.google.com/optimization) to solve Linear Programs that appear in polyhedral operations. By default we use as our solver [GUROBI](https://www.gurobi.com), however one can easily select their preferred solver by modifying the appropriate lines of the source code. The following instructions assume that you have already installed OR-Tools for C++ using the instructions [here](https://developers.google.com/optimization/install/cpp).

Other dependencies:
* CMake (latest version recommended)
* Eigen

### Install
To install the C++ library run the following commands from `/path-to-cis2m/cpp/`:
```
mkdir build
cd build
cmake ..
make
sudo make install
```

### Test
Execute the command `ctest` to make sure that the library is running as intended. 

### Usage:
After successfully installing, the library is used by the following include statement `#include "cis_generator.hpp"`.