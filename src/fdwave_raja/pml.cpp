#include "pml.hpp"

#include <math.h>

#include "constants.hpp"

// For debug
#include <stdlib.h>
#include <stdio.h>

/**
 * @param profile has dimension [i_min,i_max]
 */
void pml_profile_init(float *profile, llint i_min, llint i_max, llint n_first, llint n_last,
                      float scale)
{
    llint n = i_max-i_min+1;
    llint shift = i_min-1;

    llint first_beg = 1 + shift;
    llint first_end = n_first + shift;
    llint last_beg  = n - n_last+1 + shift;
    llint last_end  = n + shift;

    for (llint i = i_min; i <= i_max; ++i) {
        profile[i] = 0.f;
    }

    float tmp = scale / POW2(first_end-first_beg+1);
    for (llint i = 1; i <= first_end-first_beg+1; ++i) {
        profile[first_end-i+1] = POW2(i)*tmp;
    }

    for (llint i = 1; i <= last_end-last_beg+1; ++i) {
        profile[last_beg+i-1] = POW2(i)*tmp;
    } 
}

void pml_profile_extend(llint ny, llint nz,
                        float *eta, const float *etax, const float *etay, const float *etaz,
                        llint xbeg, llint xend, llint ybeg, llint yend, llint zbeg, llint zend)
{
    const llint n_ghost = 1;
    for (llint ix = xbeg-n_ghost; ix <= xend+n_ghost; ++ix) {
        for (llint iy = ybeg-n_ghost; iy <= yend+n_ghost; ++iy) {
            for (llint iz = zbeg-n_ghost; iz <= zend+n_ghost; ++iz) {
                eta[IDX3_eta0(ix,iy,iz)] = etax[ix] + etay[iy] + etaz[iz];
            }
        }
    }
}

void pml_profile_extend_all(llint ny, llint nz,
                            float *eta, const float *etax, const float *etay, const float *etaz,
                            llint xmin, llint xmax, llint ymin, llint ymax,
                            llint x1, llint x2, llint x5, llint x6,
                            llint y1, llint y2, llint y3, llint y4, llint y5, llint y6,
                            llint z1, llint z2, llint z3, llint z4, llint z5, llint z6)
{
    // Top.
    if (z1 != -1)
    pml_profile_extend(ny,nz,eta,etax,etay,etaz,xmin,xmax,ymin,ymax,z1,z2);
    // Bottom.
    if (z5 != -5)
    pml_profile_extend(ny,nz,eta,etax,etay,etaz,xmin,xmax,ymin,ymax,z5,z6);
    // Front.
    if ((y1!=-1) && (z3!=-3))
    pml_profile_extend(ny,nz,eta,etax,etay,etaz,xmin,xmax,y1,y2,z3,z4);
    // Back.
    if ((y6!=-6) && (z3!=-3))
    pml_profile_extend(ny,nz,eta,etax,etay,etaz,xmin,xmax,y5,y6,z3,z4);
    // Left.
    if ((x1!=-1) && (y3!=-3) && (z3!=-3))
    pml_profile_extend(ny,nz,eta,etax,etay,etaz,x1,x2,y3,y4,z3,z4);
    // Right.
    if ((x6!=-6) && (y3!=-3) && (z3!=-3))
    pml_profile_extend(ny,nz,eta,etax,etay,etaz,x5,x6,y3,y4,z3,z4);
}

void init_eta(llint nx, llint ny, llint nz, struct grid_t grid,
              float dt_sch,
              float *eta)
{
    for (llint i = -1; i < nx+1; ++i) {
        for (llint j = -1; j < ny+1; ++j) {
            for (llint k = -1; k < nz+1; ++k) {
                eta[IDX3_eta1(i,j,k)] = 0.f;
            }
        }
    }

    /* etax */
    float param = dt_sch * 3.f * vmax * logf(1000.f)/(2.f*grid.ndampx*grid.dx);
    float *etax =(float *) malloc(sizeof(float)*(nx+2));
    pml_profile_init(etax, 0, grid.nx+1, grid.ndampx, grid.ndampx, param);

    /* etay */
    param = dt_sch*3.f*vmax*logf(1000.f)/(2.f*grid.ndampy*grid.dy);
    float *etay =(float *) malloc(sizeof(float)*(ny+2));
    pml_profile_init(etay, 0, grid.ny+1, grid.ndampy, grid.ndampy, param);

    /* etaz */
    param = dt_sch*3.f*vmax*logf(1000.f)/(2.f*grid.ndampz*grid.dz);
    float *etaz =(float *) malloc(sizeof(float)*(nz+2));
    pml_profile_init(etaz, 0, grid.nz+1, grid.ndampz, grid.ndampz, param);

    (void)pml_profile_extend_all(ny, nz,
                eta, etax, etay, etaz,
                1, nx, 1, ny,
                grid.x1+1, grid.x2, grid.x5+1, grid.x6,
                grid.y1+1, grid.y2, grid.y3+1, grid.y4, grid.y5+1, grid.y6,
                grid.z1+1, grid.z2, grid.z3+1, grid.z4, grid.z5+1, grid.z6);

    free(etax);
    free(etay);
    free(etaz);
}
