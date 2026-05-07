#include <cmath>
#include <algorithm>

// ============================================================
// SIMPLE-like pressure-correction algorithm
//
// The momentum predictor first computes an intermediate velocity
// field using the current pressure. Since this velocity usually
// does not satisfy continuity, a pressure-correction equation is
// solved. The resulting pressure correction is then used to update
// pressure and velocity. The continuity residual measures the
// remaining mass-conservation error.
// ============================================================
/*
SIMPLE algorithm steps:

guess pressure
↓
solve momentum equation → get u*, v*
↓
u*, v* usually do not satisfy continuity
↓
solve pressure correction equation
↓
correct velocity and pressure
↓
repeat until continuity residual is small

*/

// Computes the predicted velocity field (u_star, v_star) from the momentum equations.
// Uses the current velocity and pressure fields.
// The predicted velocity is not guaranteed to satisfy continuity yet.
void momentumPredictor(FlowField2D& field,
                       const Grid2D& grid,
                       const Solver2D& solver)
{
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;
    double dt = solver.dt;

    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {
            double du_dx = d_dx_upwind(field.u, field.u, i, j, dx);
            double du_dy = d_dy_upwind(field.u, field.v, i, j, dy); 
            double dv_dx = d_dx_upwind(field.v, field.u, i, j, dx);
            double dv_dy = d_dy_upwind(field.v, field.v, i, j, dy);

            double conv_u = field.u[i][j] * du_dx + field.v[i][j] * du_dy;
            double conv_v = field.u[i][j] * dv_dx + field.v[i][j] * dv_dy;

            double diff_u = solver.nu * laplacian(field.u, i, j, dx, dy);
            double diff_v = solver.nu * laplacian(field.v, i, j, dx, dy);

            double dpdx = (field.p[i + 1][j] - field.p[i - 1][j]) / (2.0 * dx);
            double dpdy = (field.p[i][j + 1] - field.p[i][j - 1]) / (2.0 * dy);

            // u*, v* are predicted velocities because they are computed using the current pressure field. 
            // They are not yet guaranteed to satisfy incompressibility.
            // ∂u/∂t = −(u.∂u/∂x + v.∂u/∂y) + ν.∇^2u − 1/ρ(∂p/∂x) = RHS
            // Explicit Euler time step => u* = u + dt ⋅ RHS 

            field.u_star[i][j] = field.u[i][j] + dt * (-conv_u + diff_u - dpdx / solver.rho);
            field.v_star[i][j] = field.v[i][j] + dt * (-conv_v + diff_v - dpdy / solver.rho);
        }
    }

    // Copy boundaries approximately before pressure correction
    for (int j = 0; j < ny; j++)
    {
        field.u_star[0][j] = field.u[0][j];
        field.v_star[0][j] = field.v[0][j];
        field.u_star[nx - 1][j] = field.u[nx - 1][j];
        field.v_star[nx - 1][j] = field.v[nx - 1][j];
    }
    for (int i = 0; i < nx; i++)
    {
        field.u_star[i][0] = field.u[i][0];
        field.v_star[i][0] = field.v[i][0];
        field.u_star[i][ny - 1] = field.u[i][ny - 1];
        field.v_star[i][ny - 1] = field.v[i][ny - 1];
    }
}


// Solves the pressure-correction equation.
// The pressure correction p_corr is computed from the divergence of the predicted velocity.
// This correction will later be used to make the velocity field more divergence-free.
void solvePressureCorrection(FlowField2D& field,
                             const Grid2D& grid,
                             const Solver2D& solver)
{
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;
    double dt = solver.dt;

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            field.p_corr[i][j] = 0.0;

    
    // ∇²p' in central difference form
    // (p'_{i+1,j} - 2*p'_{i,j} + p'_{i-1,j})/dx² + (p'_{i,j+1} - 2*p'_{i,j} + p'_{i,j-1})/dy²    
    // Rearranging gives:
    // (p'_{i+1,j} + p'_{i-1,j})/dx² + (p'_{i,j+1} + p'_{i,j-1})/dy² - p'_{i,j} * (2/dx² + 2/dy²)

    double ae = 1.0 / (dx * dx);
    double aw = 1.0 / (dx * dx);
    double an = 1.0 / (dy * dy);
    double as = 1.0 / (dy * dy);
    double ap = ae + aw + an + as;
    
    // This gives ∇²p' = ae*p'_{i+1,j} + aw*p'_{i-1,j} + an*p'_{i,j+1} + as*p'_{i,j-1} - ap*p'_{i,j}

    for (int iter = 0; iter < solver.pressureIterations; iter++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            for (int j = 1; j < ny - 1; j++)
            {    
                // Pressure correction equation derived from continuity and momentum:   
                // ∇²p' = (ρ/dt) ∇·u* 

                double div_star = (field.u_star[i + 1][j] - field.u_star[i - 1][j]) / (2.0 * dx)
                                + (field.v_star[i][j + 1] - field.v_star[i][j - 1]) / (2.0 * dy);

                double rhs = solver.rho / dt * div_star;
                
                // p'_{i,j} = (ae*p'_{i+1,j} + aw*p'_{i-1,j} + an*p'_{i,j+1} + as*p'_{i,j-1} - rhs) / ap
                field.p_corr[i][j] = (ae * field.p_corr[i + 1][j]
                                    + aw * field.p_corr[i - 1][j]
                                    + an * field.p_corr[i][j + 1]
                                    + as * field.p_corr[i][j - 1]
                                    - rhs) / ap;
            }
        }
        applyPressureCorrectionBoundary(field, grid);
    }
}


// Corrects pressure and velocity using the pressure-correction field.
// 1. Compute ∇p'
// 2. Correct u*, v* using ∇p'
// 3. Update u, v, p with under-relaxation
void correctPressureAndVelocity(FlowField2D& field,
                                const Grid2D& grid,
                                const Solver2D& solver)
{
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;
    double dt = solver.dt;

    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {   
            // Compute the pressure correction gradients at the cell center using central differences.
            double dpdx_corr = (field.p_corr[i + 1][j] - field.p_corr[i - 1][j]) / (2.0 * dx);
            double dpdy_corr = (field.p_corr[i][j + 1] - field.p_corr[i][j - 1]) / (2.0 * dy);

            // u_correction = u* - (dt/ρ) ∂p'/∂x
            // corrected velocity is a combination of the predicted velocity and the pressure correction gradient.
            double u_corr = field.u_star[i][j] - dt / solver.rho * dpdx_corr;
            double v_corr = field.v_star[i][j] - dt / solver.rho * dpdy_corr;

            // new u = 70% corrected u + 30% old u (under-relaxation)
            field.u[i][j] = solver.alpha_u * u_corr + (1.0 - solver.alpha_u) * field.u[i][j];
            field.v[i][j] = solver.alpha_u * v_corr + (1.0 - solver.alpha_u) * field.v[i][j];
            field.p[i][j] = field.p[i][j] + solver.alpha_p * field.p_corr[i][j];
        }
    }
}


// Computes the RMS continuity residual over all interior cells.
// This measures how well the incompressibility condition div(u) = 0 is satisfied.
// check if ∇·u is close to zero to determine if we have a divergence-free velocity field, which is the ultimate goal of the pressure correction step.
double computeContinuityResidual(const FlowField2D& field,
                                 const Grid2D& grid)
{
    double sum = 0.0;
    int count = 0;

    for (int i = 1; i < grid.nx - 1; i++)
    {
        for (int j = 1; j < grid.ny - 1; j++)
        {
            double div = divergence(field, i, j, grid.dx, grid.dy);
            sum += div * div;
            count++;
        }
    }

    return std::sqrt(sum / std::max(1, count));
}
