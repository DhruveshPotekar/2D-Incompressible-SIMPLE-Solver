// ============================================================
// Boundary conditions for flat-plate boundary-layer test case
//
// left   : velocity inlet, u = U_inf, v = 0
// bottom : viscous no-slip wall, u = 0, v = 0
// top    : free-stream velocity boundary, u = U_inf, v = 0
// right  : pressure outlet, p = p_out, du/dx = dv/dx = 0
// ============================================================

void applyBoundaryConditions(FlowField2D& field,
                             const Grid2D& grid,
                             const Solver2D& solver)
{
    int nx = grid.nx;
    int ny = grid.ny;

    // Inlet
    for (int j = 0; j < ny; j++)
    {
        field.u[0][j] = solver.U_inf;
        field.v[0][j] = 0.0; 
        field.p[0][j] = field.p[1][j]; // dp/dx = 0 at inlet
    }

    // Bottom wall
    for (int i = 0; i < nx; i++)
    {
        field.u[i][0] = 0.0; // No-slip condition
        field.v[i][0] = 0.0; // No-slip condition
        field.p[i][0] = field.p[i][1]; // dp/dy = 0 at wall
    }

    // Top free stream, forces the top boundary to remain free stream 
    for (int i = 0; i < nx; i++)
    {
        field.u[i][ny - 1] = solver.U_inf;
        field.v[i][ny - 1] = 0.0;
        field.p[i][ny - 1] = field.p[i][ny - 2];
    }

    // Outlet
    for (int j = 0; j < ny; j++)
    {
        field.u[nx - 1][j] = field.u[nx - 2][j];
        field.v[nx - 1][j] = field.v[nx - 2][j];
        field.p[nx - 1][j] = solver.p_out;
    }
}

void applyPressureCorrectionBoundary(FlowField2D& field,
                                     const Grid2D& grid)
{
    int nx = grid.nx;
    int ny = grid.ny;

    // Zero normal gradient at inlet, wall, and top
    for (int j = 0; j < ny; j++)
        field.p_corr[0][j] = field.p_corr[1][j];

    for (int i = 0; i < nx; i++)
    {
        field.p_corr[i][0] = field.p_corr[i][1];
        field.p_corr[i][ny - 1] = field.p_corr[i][ny - 2];
    }

    // Fixed correction at pressure outlet
    for (int j = 0; j < ny; j++)
        field.p_corr[nx - 1][j] = 0.0;
}
