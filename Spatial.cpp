#include <cmath>
#include <algorithm>

// ============================================================
// Spatial operators
// Convective terms use first-order upwind for robustness.
// Diffusion and pressure gradients use central differences.
// ============================================================

// 1st-order Upwind scheme for convection in x-direction
double d_dx_upwind(const std::vector<std::vector<double>>& phi, //phi can be any field
                   const std::vector<std::vector<double>>& u,
                   int i, int j, double dx)
{
    if (u[i][j] >= 0.0)
        return (phi[i][j] - phi[i - 1][j]) / dx; // Backward difference for positive velocity
    else
        return (phi[i + 1][j] - phi[i][j]) / dx; // Forward difference for negative velocity
}

// 1st-order Upwind scheme for convection in y-direction
double d_dy_upwind(const std::vector<std::vector<double>>& phi,
                   const std::vector<std::vector<double>>& v,
                   int i, int j, double dy)
{
    if (v[i][j] >= 0.0)
        return (phi[i][j] - phi[i][j - 1]) / dy; // Backward difference for positive velocity
    else
        return (phi[i][j + 1] - phi[i][j]) / dy; // Forward difference for negative velocity
}

// Standard second-order central difference formula for the Laplacian operator, used in diffusion and pressure Poisson equation.
double laplacian(const std::vector<std::vector<double>>& phi,
                 int i, int j, double dx, double dy)
{
    return (phi[i + 1][j] - 2.0 * phi[i][j] + phi[i - 1][j]) / (dx * dx)
         + (phi[i][j + 1] - 2.0 * phi[i][j] + phi[i][j - 1]) / (dy * dy);
}

// Divergence of the velocity field, used in the continuity equation and pressure correction step.
double divergence(const FlowField2D& field,
                  int i, int j, double dx, double dy)
{
    double dudx = (field.u[i + 1][j] - field.u[i - 1][j]) / (2.0 * dx);
    double dvdy = (field.v[i][j + 1] - field.v[i][j - 1]) / (2.0 * dy);
    return dudx + dvdy;
}
