#include "grid.hpp"

#include "constants.hpp"
#include <stdio.h>

struct grid_t init_grid(llint nx, llint ny, llint nz, llint tsx, llint tsy,
                        llint ngpu)
{
    struct grid_t grid;
    grid.nx = nx; grid.ny = ny; grid.nz = nz;
    grid.dx = 20;  grid.dy = 20;  grid.dz = 20;
    grid.lx = 4; grid.ly = 4; grid.lz = 4;
    grid.ntaperx = 3; grid.ntapery = 3; grid.ntaperz = 3;

    const float lambdamax = vmax/_fmax;
    grid.ndampx = grid.ntaperx * lambdamax / grid.dx;
    grid.ndampy = grid.ntapery * lambdamax / grid.dy;
    grid.ndampz = grid.ntaperz * lambdamax / grid.dz;

    grid.x1 = 0;
    grid.x2 = grid.ndampx;
    grid.x3 = grid.ndampx;
    grid.x4 = grid.nx-grid.ndampx;
    grid.x5 = grid.nx-grid.ndampx;
    grid.x6 = grid.nx;

    grid.y1 = 0;
    grid.y2 = grid.ndampy;
    grid.y3 = grid.ndampy;
    grid.y4 = grid.ny-grid.ndampy;
    grid.y5 = grid.ny-grid.ndampy;
    grid.y6 = grid.ny;

    grid.z1 = 0;
    grid.z2 = grid.ndampz;
    grid.z3 = grid.ndampz;
    grid.z4 = grid.nz-grid.ndampz;
    grid.z5 = grid.nz-grid.ndampz;
    grid.z6 = grid.nz;

    grid.tsx = tsx;
    grid.tsy = tsy;
    grid.ntx = nx/tsx;
    grid.nty = ny/tsy;

    // For multi-gpu targets
    grid.ngpu = ngpu;

    printf("ndamp = %lld %lld %lld\n", grid.ndampx, grid.ndampy, grid.ndampz);
    return grid;
}
