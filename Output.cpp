#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

// ============================================================
// Tecplot and residual output
// ============================================================

void writeTecplot(const FlowField2D& field,
                  const Grid2D& grid,
                  const std::string& filename)
{
    std::ofstream file(filename);
    file << std::scientific << std::setprecision(8);

    file << "TITLE = \"2D Incompressible Flat Plate SIMPLE-like Solver\"\n";
    file << "VARIABLES = \"x\" \"y\" \"u\" \"v\" \"p\" \"velocity_mag\"\n";
    file << "ZONE I=" << grid.nx << ", J=" << grid.ny << ", F=POINT\n";

    for (int j = 0; j < grid.ny; j++)
    {
        for (int i = 0; i < grid.nx; i++)
        {
            double mag = std::sqrt(field.u[i][j] * field.u[i][j]
                                 + field.v[i][j] * field.v[i][j]);

            file << grid.x[i] << " "
                 << grid.y[j] << " "
                 << field.u[i][j] << " "
                 << field.v[i][j] << " "
                 << field.p[i][j] << " "
                 << mag << "\n";
        }
    }
}

void writeResidualHeader(const std::string& filename)
{
    std::ofstream file(filename);
    file << "TITLE = \"Residual History\"\n";
    file << "VARIABLES = \"Iteration\" \"ContinuityResidual\"\n";
    file << "ZONE F=POINT\n";
}

void appendResidual(const std::string& filename,
                    int iteration,
                    double residual)
{
    std::ofstream file(filename, std::ios::app);
    file << iteration << " " << residual << "\n";
}
