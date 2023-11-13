# bitarrays

A python interface to a C++ implementation of bitsets and bitarrays with multidimensional transforms, as in [Udovenko 2019](https://link.springer.com/chapter/10.1007/978-3-030-92062-3_12).

## compile

```sh
g++ -O3 -Wall -shared -std=c++2a -DNDEBUG -funroll-loops -fvisibility=hidden -fPIC $(python3 -m pybind11 --includes) src/bitset.cpp -o bitset$(python3-config --extension-suffix)
```
