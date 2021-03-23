#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "constants.hpp"
#include "grid.hpp"
#include "pml.hpp"
#include "data_setup.hpp"
#include "target_3d.hpp"

#include "memoryManager.hpp"
#include "RAJA/RAJA.hpp"

const float _fmax = 25.0f;
const float vmin = 1500.f;
const float vmax = 4500.f;
const float cfl = 0.8f;

void init_coef(float dx, float *coefx)
{
    float dx2 = dx*dx;
    coefx[0] = -205.f/72.f/dx2;
    coefx[1] = 8.f/5.f/dx2;
    coefx[2] = -1.f/5.f/dx2;
    coefx[3] = 8.f/315.f/dx2;
    coefx[4] = -1.f/560.f/dx2;
}

float compute_dt_sch(const float *coefx, const float *coefy, const float *coefz)
{

    float ftmp = 0.f;
    ftmp += fabsf(coefx[0]) + fabsf(coefy[0]) + fabsf(coefz[0]);
    for (uint i = 1; i < 5; i++) {
        ftmp += 2.f*fabsf(coefx[i]);
        ftmp += 2.f*fabsf(coefy[i]);
        ftmp += 2.f*fabsf(coefz[i]);
    }
    return 2.f*cfl/(sqrtf(ftmp)*vmax);
}

void gaussian_source(uint nt, float dt, float *source)
{
    const float sigma = 0.6f*_fmax;
    const float tau = 1.0f;
    const float scale = 8.0f;

    for (uint it = 1; it <= nt; ++it) {
        float t = dt*(it-1);
        source[it-1] = -2.f*scale*sigma
                     *(sigma-2.f*sigma*scale*POW2(sigma*t-tau))
                     *expf(-scale*POW2(sigma*t-tau));
    }
}

void write_io(llint nx, llint ny, llint nz,
              llint lx, llint ly, llint lz,
              const float *u, uint istep)
{
    char filename_buf[32];
    snprintf(filename_buf, sizeof(filename_buf), "snapshot.it%u.n%llu.raw", istep, nz);
    FILE *snapshot_file = fopen(filename_buf, "wb");
    for (llint k = lz; k < nz+lz; ++k) {
        for (llint j = ly; j < ny+ly; ++j) {
            for (llint i = lx; i < nx+lx; ++i) {
                fwrite(&u[IDX3_l(i,j,k)], sizeof(float),1, snapshot_file);
            }
        }
    }
    /* Clean up */
    fclose(snapshot_file);
}


int main(int argc, char *argv[])
{
    llint nx = 100, ny = 100, nz = 100;
    /* Task size (supported targets only) */
    llint tsx = 10, tsy = 10;
    /* Number of GPUs (supported targets only) */
    llint ngpus = 1;
    uint nsteps = 1000;
    uint niters = 1;
    bool disable_warm_up_iter = true;
    bool finalio = false;

    /* Process arguments */
    if( argc == 1 )
        printf( "Usage:\n ./main --grid N --nsteps M --finalio\nUsing default values, no IO.\n\n" );
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i],"--grid") == 0) {
            ++i;
            nx = strtoll(argv[i],NULL,10);
            ny = nx;
            nz = nx;
        }
        else if (strcmp(argv[i],"--tsize") == 0) {
            ++i;
            tsx = strtoll(argv[i],NULL,10);
            tsy = tsx;
            printf("tsx, tsy = %lld, %lld\n", tsx, tsy);
        }
        else if (strcmp(argv[i],"--nsteps") == 0) {
            ++i;
            nsteps = strtoll(argv[i],NULL,10);
            printf("nsteps = %d\n", nsteps);
        }
        else if (strcmp(argv[i],"--niters") == 0) {
            ++i;
            niters = strtoll(argv[i],NULL,10);
            printf("niters = %d\n", niters);
        }
        else if (strcmp(argv[i],"--warm-up") == 0) {
            disable_warm_up_iter = false;
            printf("warm up iteration is enabled\n");
        }
        else if (strcmp(argv[i],"--finalio") == 0) {
            finalio = true;
            printf("writing final wavefile is enabled\n");
        }
        else if (strcmp(argv[i],"--ngpus") == 0) {
            ++i;
            ngpus = strtoll(argv[i],NULL,10);
            printf("ngpus = %lld\n", ngpus);
        }
    }

    double time_kernel = 0.0;
    double time_modeling = 0.0;
    bool warm_up_iter = !disable_warm_up_iter;



    for (uint iiter = 0; iiter < (disable_warm_up_iter ? niters : niters+1); iiter++) {

        const struct grid_t grid = init_grid(nx,ny,nz,tsx,tsy,ngpus);

        printf("grid = %lld %lld %lld\n", grid.nx, grid.ny, grid.nz);

        const llint lx = grid.lx, ly = grid.ly, lz = grid.lz;

        const llint sx = nx/2, sy = ny/2, sz = nz/2;

        //
        // Define range indices
        //
        RAJA::TypedRangeSegment<llint> KLRange(-lx, nx+lx-1);
        RAJA::TypedRangeSegment<llint> JLRange(-ly, ny+ly-1);
        RAJA::TypedRangeSegment<llint> ILRange(-lz, nz+lz-1);
        //
        RAJA::TypedRangeSegment<llint> XvpRange(0, nx-1);
        RAJA::TypedRangeSegment<llint> YvpRange(0, ny-1);
        RAJA::TypedRangeSegment<llint> ZvpRange(0, nz-1);

        RAJA::TypedRangeSegment<llint> XphiRange(0, nx-1);
        RAJA::TypedRangeSegment<llint> YphiRange(0, ny-1);
        RAJA::TypedRangeSegment<llint> ZphiRange(0, nz-1);

        RAJA::TypedRangeSegment<llint> XetaRange(-1, nx);
        RAJA::TypedRangeSegment<llint> YetaRange(-1, ny);
        RAJA::TypedRangeSegment<llint> ZetaRange(-1, nz);

        //
        // Define the nested loop
        //
        using KJI_EXECPOL = RAJA::KernelPolicy<
                            RAJA::statement::For<2, RAJA::loop_exec,    // k
                            RAJA::statement::For<1, RAJA::loop_exec,  // j
                            RAJA::statement::For<0, RAJA::loop_exec,// i 
                            RAJA::statement::Lambda<0>
                            > 
                          > 
                        > 
                      >;

        //float *u = (float *) malloc(sizeof(float)*(nx+2*lx)*(ny+2*ly)*(nz+2*lz));
        //float *v = (float *)malloc(sizeof(float)*(nx+2*lx)*(ny+2*ly)*(nz+2*lz));
        float *u = memoryManager::allocate<float>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
        float *v = memoryManager::allocate<float>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
        

        // PML arrays
        //float *phi =(float *) malloc(sizeof(float)*nx*ny*nz);
        //float *eta =(float *) malloc(sizeof(float)*(nx+2)*(ny+2)*(nz+2));
        float *phi = memoryManager::allocate<float>(nx*ny*nz);
        float *eta = memoryManager::allocate<float>((nx+2)*(ny+2)*(nz+2)); 

        // Wave field initialization
        //for (llint i = -lx; i < nx+lx; ++i) {
        //    for (llint j = -ly; j < ny+ly; ++j) {
        //        for (llint k = -lz; k < nz+lz; ++k) {
        //            u[IDX3_l(i,j,k)] = 0.0f;
        //            v[IDX3_l(i,j,k)] = 0.0f;
        //        }
        //    }
        //}
        RAJA::kernel<KJI_EXECPOL>( RAJA::make_tuple(ILRange, JLRange, KLRange),
                                   [=] (llint i, llint j, llint k) 
                                   { 
                                      //printf( " (%d, %d, %d) \n", (int)(i), (int)(j), (int)(k));
                                      u[IDX3_l(i,j,k)] = 0.0f;
                                      v[IDX3_l(i,j,k)] = 0.0f;
                                   }
                                 );

        //float *source =(float *) malloc(sizeof(float)*nsteps);
        //float *vp =(float *) malloc(sizeof(float)*nx*ny*nz);
        float *source = memoryManager::allocate<float>(nsteps);
        float *vp     = memoryManager::allocate<float>(nx*ny*nz);

        float coefx[5], coefy[5], coefz[5];
        init_coef(grid.dx, coefx);
        init_coef(grid.dy, coefy);
        init_coef(grid.dz, coefz);

        const float dt_sch = compute_dt_sch(coefx, coefy, coefz);
        // Init vp and phi
        const float vp_all = POW2(2000.f*dt_sch);
        //for (llint i = 0; i < nx; ++i) {
        //    for (llint j = 0; j < ny; ++j) {
        //        for (llint k = 0; k < nz; ++k) {
        //            phi[IDX3(i,j,k)] = 0.f;
        //            vp[IDX3(i,j,k)] = vp_all;
        //        }
        //    }
        //}
        RAJA::kernel<KJI_EXECPOL>( RAJA::make_tuple(XvpRange, YvpRange, ZvpRange),
                                   [=] (llint i, llint j, llint k) 
                                   { 
                                    //phi[IDX3(i,j,k)] = 0.f;
                                    vp[IDX3(i,j,k)] = vp_all;
                                   }
                                 );

        RAJA::kernel<KJI_EXECPOL>( RAJA::make_tuple(XvpRange, YvpRange, ZvpRange),
                                   [=] (llint i, llint j, llint k) 
                                   { 
                                    //phi[IDX3(i,j,k)] = 0.f;
                                    phi[IDX3(i,j,k)] = 0;
                                   }
                                 );        
        (void) gaussian_source(nsteps,dt_sch,source);
        // Init PML
        init_eta(nx, ny, nz, grid, dt_sch,
                 eta);

        const float hdx_2 = 1.f / (4.f * POW2(grid.dx));
        const float hdy_2 = 1.f / (4.f * POW2(grid.dy));
        const float hdz_2 = 1.f / (4.f * POW2(grid.dz));

        struct timespec start,end,start_m,end_m;
        const uint npo = 100;
        const bool l_snapshot = false;
        const uint nsnapshot_freq = 100;

        //target_init(grid,nsteps,u,v,phi,eta,coefx,coefy,coefz,vp,source);

        // Time loop
        clock_gettime(CLOCK_REALTIME, &start_m);
        // #ifdef OMP_PARALLEL_MASTER_IN_MAIN
        //     #pragma omp parallel default(shared)
        //     #pragma omp master
        // #endif
        #ifdef __NVCC__
        target(
            nsteps, &time_kernel,
            nx, ny, nz,
            grid.x1, grid.x2, grid.x3, grid.x4, grid.x5, grid.x6,
            grid.y1, grid.y2, grid.y3, grid.y4, grid.y5, grid.y6,
            grid.z1, grid.z2, grid.z3, grid.z4, grid.z5, grid.z6,
            lx, ly, lz,
            sx, sy, sz,
            hdx_2, hdy_2, hdz_2,
            coefx, coefy, coefz,
            u, v, vp,
            phi, eta, source
        );
        if (warm_up_iter) {
            time_kernel = 0;
        }
        #else
        for (uint istep = 1; istep <= nsteps; ++istep) {

            /* Write snapshot if requested */
            if (l_snapshot && istep % nsnapshot_freq == 0) {
                write_io(nx, ny, nz, lx, ly, lz,
                         u, istep);
            }

            /*******************************
             * Main loop
             *******************************/
            clock_gettime(CLOCK_REALTIME, &start);
            target_3d( nx, ny, nz,
                       grid.x1, grid.x2, grid.x3, grid.x4, grid.x5, grid.x6,
                       grid.y1, grid.y2, grid.y3, grid.y4, grid.y5, grid.y6,
                       grid.z1, grid.z2, grid.z3, grid.z4, grid.z5, grid.z6,
                       lx, ly, lz,
                       hdx_2, hdy_2, hdz_2,
                       coefx, coefy, coefz,
                       u, v, vp, phi, eta );
            clock_gettime(CLOCK_REALTIME, &end);
            if (!warm_up_iter) {
                time_kernel += (end.tv_sec  - start.tv_sec) +
                               (double)(end.tv_nsec - start.tv_nsec) / 1.0e9;
            }

            // Print out
            if (istep % npo == 0) {
                printf("time step %u / %u\n", istep, nsteps);

                #ifdef MINMAX_U_IN_TIME_LOOP
                    float min_v, max_v;
                    find_min_max_u(grid, v, &min_v, &max_v);
                    printf("min_v,  max_v = %f, %f\n", min_v, max_v);
                #endif
            }

            // Add source
            kernel_add_source(grid, v, source, istep, sx, sy, sz);

            // Swap u and v
            float *t = u;
            u = v;
            v = t;
        }
        #endif
        clock_gettime(CLOCK_REALTIME, &end_m);
        if (!warm_up_iter) {
            time_modeling += (end_m.tv_sec  - start_m.tv_sec) +
                             (double)(end_m.tv_nsec - start_m.tv_nsec) / 1.0e9;
        }

        float min_u, max_u;
        find_min_max_u(grid, u, &min_u, &max_u);
        printf("FINAL min_u,  max_u = %f, %f\n", min_u, max_u);

        if( finalio )
            write_io(nx, ny, nz, lx, ly, lz, u, nsteps);

        //target_finalize(grid,nsteps,u,v,phi,eta,coefx,coefy,coefz,vp,source);
        free(u);
        free(v);
        free(phi);
        free(eta);
        free(source);
        free(vp);

        warm_up_iter = false;
    }

    printf("Time kernel: %g s\n", time_kernel / niters);
    printf("Time modeling: %g s\n", time_modeling / niters);
}
