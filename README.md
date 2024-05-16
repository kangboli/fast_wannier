# A Fortran Implementation of the Truncated Density Convolution

## Build and Run

Prereq: `OpenBLAS` or `MKL`.

```sh
mkdir build
cd build
cmake ..
make

./wannier.x ../data/Si
```

This will take many iterations because it is running a fixed step gradient
descent. There is a number of other examples problems in the `data` directory.

##  Files

- The main code is in `oracles.f90`, which computes the total spread and its gradient.
- `oracles.f90` uses `zgpadm.f90`, which is downloaded from
  <https://www.maths.uq.edu.au/expokit/>. This is the matrix exponential
  through the Padé approximation.
- `optimizer.f90`, `param_parser.f90`, and `main.f90` loads and run some small
  problems for the purpose of showing how to use `oracles.f90`.

