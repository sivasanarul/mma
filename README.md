# MMA C++

A compact C++17 implementation of the Method of Moving Asymptotes (MMA) for
constrained nonlinear optimization. The repository includes the canonical
three-variable toy problem, a dense partial-pivoting linear solver, CMake build
support, and regression tests.

The current code implements MMA. It does **not** yet implement the complete
globally convergent GCMMA workflow (`gcmmasub`, `concheck`, and `raaupdate`).
The current dense Newton-system path requires the number of constraints to be
smaller than the number of design variables (`0 < m < n`).

## Provenance

The optimizer is a C++ port of the GPL-licensed MATLAB MMA/GCMMA distribution
written by Krister Svanberg. Its interface, variable names, primal-dual Newton
flow, and toy problem also correspond to the GPL Python port maintained by
Arjen Deetman. The exact source used for the original 2019 C++ port is not
recorded in the historical commits, so both parts of the public code lineage
are credited in [NOTICE](NOTICE).

The previous `LU.cpp` implementation matched third-party example code whose
license was not stated. It has been removed and replaced with a newly written
partial-pivoting Gaussian-elimination solver.

## Build and run

Using CMake 3.16 or newer:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/mma_demo
```

Or compile directly with a C++17 compiler:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -I. \
  main.cpp mma.cpp toy2.cpp linear_solver.cpp -o mma_demo
./mma_demo
```

Expected convergence summary:

```text
MMA done after 16 iterations
KKT norm: 4.28402e-08
2.01752  1.78001  1.23751
```

The final digits can vary slightly across compilers and platforms.

## Tests

With CMake:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The tests check the pivoted dense solver and require the toy MMA problem to
reach its KKT convergence tolerance.

## Correctness repairs

The modernization fixes several defects in the original port:

- the outer iteration counter now advances;
- the calculated KKT norm controls termination;
- KKT residual accumulators are initialized;
- Newton work arrays use problem-dependent dimensions;
- the dual `dzet` update includes the missing `zet` multiplier;
- line-search trial points are formed from the saved iterate; and
- the unlicensed, non-pivoting Crout routine is replaced by a tested solver.

## References

- K. Svanberg, “The Method of Moving Asymptotes — A New Method for
  Structural Optimization,” *International Journal for Numerical Methods in
  Engineering*, 24(2), 359–373, 1987.
  <https://doi.org/10.1002/nme.1620240207>
- K. Svanberg, “MMA and GCMMA — Two Methods for Nonlinear Optimization,”
  2007. <https://people.kth.se/~krille/gcmma07.pdf>
- Original GPL MATLAB distribution: <https://www.smoptit.se/>
- A. Deetman, *GCMMA-MMA-Python*:
  <https://github.com/arjendeetman/GCMMA-MMA-Python>

## License

This derivative work is distributed under the GNU General Public License,
version 3 or any later version. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
