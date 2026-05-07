#include <vector>

// ============================================================
// 2D structured Cartesian grid, cell-centred layout
// ============================================================

//Struct = data

//Function = logic that builds/manipulates that data


// A struct is like a container that holds related data together
// It does not execute anything. It only describes what a grid is.
struct Grid2D
{
    int nx, ny;
    double x_min, x_max;
    double y_min, y_max;
    double dx, dy;
    std::vector<double> x;  // 1D dynamic array to store x-coordinates of cell centers
    std::vector<double> y;  // 1D dynamic array to store y-coordinates of cell centers
};

// This is a function that creates and returns a Grid2D object. 
// It takes the number of grid points in x and y directions, and the physical dimensions of the domain as input parameters. 
// It calculates the grid spacing (dx and dy) and fills the x and y coordinate arrays with the positions of the cell centers.
Grid2D createGrid(int nx, int ny,
                  double x_min, double x_max,
                  double y_min, double y_max)
{
    Grid2D grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.x_min = x_min;
    grid.x_max = x_max;
    grid.y_min = y_min;
    grid.y_max = y_max;

    grid.dx = (x_max - x_min) / (nx - 1);
    grid.dy = (y_max - y_min) / (ny - 1);

    grid.x.resize(nx); // Resize the x vector to hold nx values
    grid.y.resize(ny); // Resize the y vector to hold ny values

    for (int i = 0; i < nx; i++)
        grid.x[i] = x_min + i * grid.dx;

    for (int j = 0; j < ny; j++)
        grid.y[j] = y_min + j * grid.dy;

    return grid;
}
