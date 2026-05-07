#include <iostream>
#include <chrono>
#include <string>

#include "Grid.cpp"
#include "Field.cpp"
#include "Solver.cpp"
#include "Initial.cpp"
#include "Spatial.cpp"
#include "Boundary.cpp"
#include "SIMPLE.cpp"
#include "Output.cpp"

int main()
{
    // ============================================================
    // Flat-plate domain
    // x: streamwise direction, y: wall-normal direction
    // ============================================================

    int nx = 161;
    int ny = 81;

    double x_min = 0.0;
    double x_max = 0.50;   // 0.5 m plate length
    double y_min = 0.0;
    double y_max = 0.05;   // 5 cm domain height

    Grid2D grid = createGrid(nx, ny, x_min, x_max, y_min, y_max);
    FlowField2D field = createFlowField(nx, ny);
    Solver2D solver = createSolver();

    solver.dt = computeTimeStep(grid, solver);

    initializeFlow(field, grid, solver);
    applyBoundaryConditions(field, grid, solver);

    std::string outputDir = "results/PI_100/";
    writeResidualHeader(outputDir + "residual.dat");
    writeTecplot(field, grid, outputDir + "solution_initial.dat");

    std::cout << "2D incompressible viscous flat-plate solver\n";
    std::cout << "U_inf = " << solver.U_inf << " m/s\n";
    std::cout << "rho   = " << solver.rho << " kg/m3\n";
    std::cout << "mu    = " << solver.mu << " Pa.s\n";
    std::cout << "dt    = " << solver.dt << " s\n";

    auto start = std::chrono::high_resolution_clock::now();

    double residual = 1.0;
    int iter = 0;

    while (iter < solver.maxIter && residual > solver.tolerance)
    {
        applyBoundaryConditions(field, grid, solver);

        momentumPredictor(field, grid, solver);
        solvePressureCorrection(field, grid, solver);
        correctPressureAndVelocity(field, grid, solver);

        applyBoundaryConditions(field, grid, solver);

        residual = computeContinuityResidual(field, grid);
        appendResidual(outputDir + "residual.dat", iter, residual);

        if (iter % 500 == 0)
        {
            std::cout << "Iteration: " << iter
                      << "  Continuity residual: " << residual << "\n";
        }

        if (iter % 5000 == 0 && iter > 0)
        {
            std::string name = outputDir + "solution_" + std::to_string(iter) + ".dat";
            writeTecplot(field, grid, name);
        }

        iter++;
    }

    writeTecplot(field, grid, outputDir + "solution_final.dat");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Finished.\n";
    std::cout << "Iterations = " << iter << "\n";
    std::cout << "Final residual = " << residual << "\n";
    std::cout << "Wall time = " << elapsed.count() << " s\n";

    return 0;
}
