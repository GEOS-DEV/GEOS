#include "data_setup.hpp"

#include <float.h>
#include <math.h>

//void target_init(struct grid_t grid, uint nsteps,
//                 const float *u, const float *v, const float *phi,
//                 const float *eta, const float *coefx, const float *coefy,
//                 const float *coefz, const float *vp, const float *source)
//{
//    // Nothing needed
//}

//void target_finalize(struct grid_t grid, uint nsteps,
//                     const float *u, const float *v, const float *phi,
//		     const float *eta, const float *coefx, const float *coefy,
//                     const float *coefz, const float *vp, const float *source)
//{
//    // Nothing needed
//}

void kernel_add_source(struct grid_t grid,
                       float *u, const float *source, llint istep,
                       llint sx, llint sy, llint sz)
{
    const llint ny = grid.ny;
    const llint nz = grid.nz;
    const llint lx = grid.lx;
    const llint ly = grid.ly;
    const llint lz = grid.lz;
    u[IDX3_l(sx,sy,sz)] += source[istep-1];
}

void find_min_max_u(struct grid_t grid,
                    const float *u, float *min_u, float *max_u)
{
    const llint nx = grid.nx;
    const llint ny = grid.ny;
    const llint nz = grid.nz;
    const llint lx = grid.lx;
    const llint ly = grid.ly;
    const llint lz = grid.lz;

    *min_u = FLT_MAX;
    *max_u = FLT_MIN;
    for (llint i = -lx; i < nx+lx; ++i) {
        for (llint j = -ly; j < ny+ly; ++j) {
            for (llint k = -lz; k < nz+lz; ++k) {
                *min_u = fminf(*min_u, u[IDX3_l(i,j,k)]);
                *max_u = fmaxf(*max_u, u[IDX3_l(i,j,k)]);
            }
        }
    }
}
