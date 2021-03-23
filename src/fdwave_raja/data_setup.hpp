#ifndef DATA_SETUP_HPP
#define DATA_SETUP_HPP

#include "constants.hpp"
#include "grid.hpp"

//void target_init(struct grid_t grid, uint nsteps,
//                 const float *u, const float *v, const float *phi,
//                 const float *eta, const float *coefx, const float *coefy,
//                 const float *coefz, const float *vp, const float *source);

//void target_finalize(struct grid_t grid, uint nsteps,
//                     const float *u, const float *v, const float *phi,
//                     const float *eta, const float *coefx, const float *coefy,
//                     const float *coefz, const float *vp, const float *source);

void kernel_add_source(struct grid_t grid,
                       float *u, const float *source, llint istep,
                       llint sx, llint sy, llint sz);

void find_min_max_u(struct grid_t grid,
                    const float *u, float *min_u, float *max_u);

#endif
