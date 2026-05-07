# 2D Incompressible SIMPLE Solver

This repository contains a 2D incompressible Navier–Stokes solver written in C++ for a flat-plate boundary-layer test case. The solver is based on a pressure-correction (SIMPLE-like) approach on a structured Cartesian grid.

The code is developed for learning and understanding the fundamentals of pressure-based CFD solvers.

---

## Governing Equations

The solver advances the incompressible laminar Navier–Stokes equations:

- Continuity equation:
  
  ∇ · u = 0

- Momentum equations:
  
  ∂u/∂t + (u · ∇)u = - (1/ρ) ∇p + ν ∇²u

where:
- ρ = density  
- ν = kinematic viscosity  

---

## Numerical Method

- Finite-difference discretization on a structured grid
- First-order upwind scheme for convective terms
- Second-order central differencing for diffusion and pressure gradients
- Explicit time stepping (pseudo-transient approach)
- SIMPLE-like pressure correction method:
  - Momentum predictor → u*, v*
  - Pressure correction (Poisson equation)
  - Velocity and pressure update with under-relaxation

---

## Boundary Conditions

- **Inlet (left):**  
  u = 10 m/s, v = 0  

- **Wall (bottom):**  
  No-slip condition: u = 0, v = 0  

- **Top boundary:**  
  Free-stream velocity: u = 10 m/s, v = 0  

- **Outlet (right):**  
  Fixed gauge pressure: p = 0  
  Zero-gradient velocity  

---

## Compilation

```bash
g++ main.cpp -O2 -std=c++17 -o simple_flat_plate

---

## Convergence Study

A comparison was performed for different numbers of pressure-correction iterations
(100, 200, 300) per SIMPLE step.

- Increasing the number of pressure iterations improves the accuracy of the pressure-correction Poisson solve.
- This leads to smoother convergence and fewer oscillations in the continuity residual.
- However, a higher number of pressure iterations also increases the computational cost per SIMPLE iteration.
- Therefore, an optimal value should balance convergence stability and total runtime.
- All cases converge toward the same final residual level, confirming numerical consistency.

