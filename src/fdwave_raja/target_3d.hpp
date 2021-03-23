#ifndef TARGET_3D_HPP
#define TARGET_3D_HPP

#include "constants.hpp"
#include "grid.hpp"

void target_3d(llint nx, llint ny, llint nz,
               llint x1, llint x2, llint x3, llint x4, llint x5, llint x6,
               llint y1, llint y2, llint y3, llint y4, llint y5, llint y6,
               llint z1, llint z2, llint z3, llint z4, llint z5, llint z6,
               llint lx, llint ly, llint lz,
               float hdx_2, float hdy_2, float hdz_2,
               const float *__restrict__ coefx, const float *coefy, const float *coefz,
               const float *u, float *v, const float *vp,
               float *phi, const float *eta);

#endif
