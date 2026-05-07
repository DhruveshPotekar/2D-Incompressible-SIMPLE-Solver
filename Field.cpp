#include <vector>

// ============================================================
// Flow-field storage for 2D incompressible flow
// ============================================================

// This defines what a flow field consists of.
// It stores all primary and auxiliary variables, each variable is a 2D field over the grid.
struct FlowField2D
{
    std::vector<std::vector<double>> u;       // x-velocity
    std::vector<std::vector<double>> v;       // y-velocity
    std::vector<std::vector<double>> p;       // pressure, gauge pressure

    std::vector<std::vector<double>> u_star;  // predicted x-velocity
    std::vector<std::vector<double>> v_star;  // predicted y-velocity
    std::vector<std::vector<double>> p_corr;  // pressure correction
};

// It is a helper function, creates a 2D array of size nx × ny, initializes all entries to zero.
std::vector<std::vector<double>> createArray(int nx, int ny, double value = 0.0) 
{
    return std::vector<std::vector<double>>(nx, std::vector<double>(ny, value));
}

// It is a factory function for your flow field. 
// It initializes all the fields (u, v, p, u_star, v_star, p_corr) to zero and returns a FlowField2D object.
FlowField2D createFlowField(int nx, int ny)
{
    FlowField2D field;
    field.u      = createArray(nx, ny, 0.0);
    field.v      = createArray(nx, ny, 0.0);
    field.p      = createArray(nx, ny, 0.0);
    field.u_star = createArray(nx, ny, 0.0);
    field.v_star = createArray(nx, ny, 0.0);
    field.p_corr = createArray(nx, ny, 0.0);
    return field;
}
