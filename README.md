# 2D Incompressible Flat-Plate Solver

This is an educational C++ solver for a 2D incompressible viscous flat-plate boundary-layer case. It follows the same modular style as the earlier 2D advection solver.

## Physics

The solver advances the incompressible laminar Navier-Stokes equations with constant density and viscosity:

- continuity equation
- x-momentum equation
- y-momentum equation

A SIMPLE-like pressure-correction loop is used to couple pressure and velocity.

## Boundary conditions

- inlet: `u = 10 m/s`, `v = 0`
- bottom wall: no slip, `u = 0`, `v = 0`
- top boundary: free-stream velocity, `u = 10 m/s`, `v = 0`
- outlet: fixed gauge pressure, `p = 0`, with zero-gradient velocity

## Compile

```bash
g++ main.cpp -O2 -std=c++17 -o simple_flat_plate
```

## Run

```bash
./simple_flat_plate
```

## Output

The code writes Tecplot ASCII files:

- `solution_initial.dat`
- `solution_final.dat`
- intermediate solution files every 1000 iterations
- `residual.dat`

## Note

This is not yet a production-quality CFD solver. It is a first learning implementation. For a stronger finite-volume SIMPLE solver, the next step is to move to staggered or Rhie-Chow stabilized collocated arrangement.
