// ============================================================
// Initial condition: uniform free stream with no-slip wall applied later
// ============================================================

// This function initializes the flow field with a uniform free-stream velocity and zero pressure.
void initializeFlow(FlowField2D& field,
                    const Grid2D& grid,
                    const Solver2D& solver)
{
    for (int i = 0; i < grid.nx; i++)
    {
        for (int j = 0; j < grid.ny; j++)
        {
            field.u[i][j] = solver.U_inf;
            field.v[i][j] = 0.0;
            field.p[i][j] = solver.p_out; // Gauge pressure initialized to zero
        }
    }
}
