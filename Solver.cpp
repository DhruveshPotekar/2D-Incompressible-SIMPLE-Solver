#include <cmath>
#include <algorithm>

// ============================================================
// Solver and physical parameters
// ============================================================

// This struct defines the parameters for the solver. 
// Just a blueprint, it does not execute anything. It only describes what a solver is.
struct Solver2D
{
    double rho;       // density
    double mu;        // dynamic viscosity
    double nu;        // kinematic viscosity

    double U_inf;     // inlet / free-stream velocity
    double p_out;     // outlet gauge pressure

    double cfl;
    double dt;

    double alpha_u;   // velocity under-relaxation
    double alpha_p;   // pressure under-relaxation

    int maxIter;
    int pressureIterations;
    double tolerance;
};

// A function that creates and returns a Solver2D object
Solver2D createSolver()
{
    Solver2D solver;

    solver.rho = 1.225;
    solver.mu  = 1.789e-5;
    solver.nu  = solver.mu / solver.rho;

    solver.U_inf = 10.0;
    solver.p_out = 0.0;        

    solver.cfl = 0.9;
    solver.dt = 0.0;

    solver.alpha_u = 0.7;
    solver.alpha_p = 0.3;

    solver.maxIter = 20000; //10000
    solver.pressureIterations = 100; // should try 200 and 300 for better convergence, but it will increase runtime significantly
    solver.tolerance = 1.0e-4;

    return solver;
}

// A stability constraint function, as the solver is explicit in time, 
// we need to compute a stable time step size based on the grid spacing and flow parameters.
double computeTimeStep(const Grid2D& grid, const Solver2D& solver)
{   
    // Compute the time step based on the CFL condition for convection and diffusion stability.
    double dt_conv = solver.cfl * std::min(grid.dx, grid.dy) / solver.U_inf;
    double dt_diff = 0.25 * std::min(grid.dx * grid.dx, grid.dy * grid.dy) / solver.nu;
    return std::min(dt_conv, dt_diff);
}dfd
