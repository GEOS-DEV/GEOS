#include "target_3d.hpp"
#include "grid.hpp"
#include "constants.hpp"

#include "memoryManager.hpp"
#include "RAJA/RAJA.hpp"

#define LAP (coef0*u[IDX3_l(i,j,k)] \
                +coefx[1]*(u[IDX3_l(i+1,j,k)]+u[IDX3_l(i-1,j,k)]) \
                +coefy[1]*(u[IDX3_l(i,j+1,k)]+u[IDX3_l(i,j-1,k)]) \
                +coefz[1]*(u[IDX3_l(i,j,k+1)]+u[IDX3_l(i,j,k-1)]) \
                +coefx[2]*(u[IDX3_l(i+2,j,k)]+u[IDX3_l(i-2,j,k)]) \
                +coefy[2]*(u[IDX3_l(i,j+2,k)]+u[IDX3_l(i,j-2,k)]) \
                +coefz[2]*(u[IDX3_l(i,j,k+2)]+u[IDX3_l(i,j,k-2)]) \
                +coefx[3]*(u[IDX3_l(i+3,j,k)]+u[IDX3_l(i-3,j,k)]) \
                +coefy[3]*(u[IDX3_l(i,j+3,k)]+u[IDX3_l(i,j-3,k)]) \
                +coefz[3]*(u[IDX3_l(i,j,k+3)]+u[IDX3_l(i,j,k-3)]) \
                +coefx[4]*(u[IDX3_l(i+4,j,k)]+u[IDX3_l(i-4,j,k)]) \
                +coefy[4]*(u[IDX3_l(i,j+4,k)]+u[IDX3_l(i,j-4,k)]) \
                +coefz[4]*(u[IDX3_l(i,j,k+4)]+u[IDX3_l(i,j,k-4)]))

void compute_inner_3d(RAJA::Index_type const x3, const RAJA::Index_type x4, const RAJA::Index_type y3,
                      RAJA::Index_type const y4, const RAJA::Index_type z3, const RAJA::Index_type z4,
                      RAJA::Index_type const lx, const RAJA::Index_type ly, const RAJA::Index_type lz,
                      RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefxViewConst,
                      RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefyViewConst,
                      RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefzViewConst,                    
                      RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ uViewConst,
                      RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ vpViewConst,
                      RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ vView )
{
  const float coef0 = coefxViewConst(0) + coefyViewConst(0) + coefzViewConst(0);

  RAJA::RangeSegment XRange(x3, x4);
  RAJA::RangeSegment YRange(y3, y4);
  RAJA::RangeSegment ZRange(z3, z4);
  
  using KJI_EXECPOL = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::loop_exec,
    RAJA::statement::For<1, RAJA::loop_exec,
    RAJA::statement::For<2, RAJA::loop_exec,
    RAJA::statement::Lambda<0> > > > >;

    RAJA::kernel<KJI_EXECPOL>( RAJA::make_tuple(XRange, YRange, ZRange),
                               [&coef0,lx,ly,lz,
                                coefxViewConst,coefyViewConst,coefzViewConst,
                                uViewConst,vpViewConst,vView] (RAJA::Index_type const i, RAJA::Index_type const j, RAJA::Index_type const k) 
                               {
                                 const float lap = (coef0*uViewConst(i+lx,j+ly,k+lz) 
                                                    +coefxViewConst(1)*(uViewConst(i+1+lx,j+ly,k+lz)+uViewConst(i-1+lx,j+ly,k+lz))
                                                    +coefyViewConst(1)*(uViewConst(i+lx,j+1+ly,k+lz)+uViewConst(i+lx,j-1+ly,k+lz))
                                                    +coefzViewConst(1)*(uViewConst(i+lx,j+ly,k+1+lz)+uViewConst(i+lx,j+ly,k-1+lz))
                                                    +coefxViewConst(2)*(uViewConst(i+2+lx,j+ly,k+lz)+uViewConst(i-2+lx,j+ly,k+lz))
                                                    +coefyViewConst(2)*(uViewConst(i+lx,j+2+ly,k+lz)+uViewConst(i+lx,j-2+ly,k+lz))
                                                    +coefzViewConst(2)*(uViewConst(i+lx,j+ly,k+2+lz)+uViewConst(i+lx,j+ly,k-2+lz))
                                                    +coefxViewConst(3)*(uViewConst(i+3+lx,j+ly,k+lz)+uViewConst(i-3+lx,j+ly,k+lz))
                                                    +coefyViewConst(3)*(uViewConst(i+lx,j+3+ly,k+lz)+uViewConst(i+lx,j-3+ly,k+lz))
                                                    +coefzViewConst(3)*(uViewConst(i+lx,j+ly,k+3+lz)+uViewConst(i+lx,j+ly,k-3+lz))
                                                    +coefxViewConst(4)*(uViewConst(i+4+lx,j+ly,k+lz)+uViewConst(i-4+lx,j+ly,k+lz))
                                                    +coefyViewConst(4)*(uViewConst(i+lx,j+4+ly,k+lz)+uViewConst(i+lx,j-4+ly,k+lz))
                                                    +coefzViewConst(4)*(uViewConst(i+lx,j+ly,k+4+lz)+uViewConst(i+lx,j+ly,k-4+lz)));
                                 vView( i+lx,j+ly,k+lz ) =
                                   2.f*uViewConst( i+lx,j+ly,k+lz )
                                   -vView( i+lx,j+ly,k+lz )
                                   +vpViewConst( i,j,k )*lap;
                               } );

}

void target_inner_3d(llint nx, llint ny, llint nz,
                     llint x3, llint x4, llint y3,
                     llint y4, llint z3, llint z4,
                     llint lx, llint ly, llint lz,
                     const float *__restrict__ coefx,
                     const float *__restrict__ coefy,
                     const float *__restrict__ coefz,
                     const float *__restrict__ u,
                     float *__restrict__ v,
                     const float *__restrict__ vp)
{
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxViewConst( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyViewConst( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzViewConst( coefz, 5 );
 
  RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > uViewConst( u, RAJA::Layout<3>(nx+2*lx, ny+2*ly, nz+2*lz) );
  RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > vView( v, RAJA::Layout<3>(nx+2*lx, ny+2*ly, nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > vpViewConst( vp, RAJA::Layout<3>(nx, ny, nz) );    
    
  compute_inner_3d(RAJA::Index_type(x3), RAJA::Index_type(x4),
                   RAJA::Index_type(y3), RAJA::Index_type(y4),
                   RAJA::Index_type(z3), RAJA::Index_type(z4),
                   RAJA::Index_type(lx), RAJA::Index_type(ly), RAJA::Index_type(lz),
                   coefxViewConst, coefyViewConst, coefzViewConst,
                   uViewConst, vpViewConst, vView );

  // RAJA_UNUSED_VAR( nx );
  // float coef0 = coefx[0] + coefy[0] + coefz[0];
  // for (llint i = x3; i < x4; ++i) {
  //   for (llint j = y3; j < y4; ++j) {
  //     for (llint k = z3; k < z4; ++k) {
  //    float lap = LAP;
  //    v[IDX3_l(i,j,k)] = 2.f*u[IDX3_l(i,j,k)]-v[IDX3_l(i,j,k)]+vp[IDX3(i,j,k)]*lap;
  //     }
  //   }
  // }
}

void compute_pml_3d(RAJA::Index_type const x3, const RAJA::Index_type x4, const RAJA::Index_type y3,
                    RAJA::Index_type const y4, const RAJA::Index_type z3, const RAJA::Index_type z4,
                    RAJA::Index_type const lx, const RAJA::Index_type ly, const RAJA::Index_type lz,
                    float const & hdx_2, float const & hdy_2, float const & hdz_2,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefxViewConst,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefyViewConst,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefzViewConst,                      
                    RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ uViewConst,
                    RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ vpViewConst,
                    RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ etaViewConst,             
                    RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ vView,
                    RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ phiView )
{
  const float coef0 = coefxViewConst(0) + coefyViewConst(0) + coefzViewConst(0);
#if 0
  RAJA::RangeSegment XRange(x3, x4);
  RAJA::RangeSegment YRange(y3, y4);
  RAJA::RangeSegment ZRange(z3, z4);
  
  using KJI_EXECPOL = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::loop_exec, 
    RAJA::statement::For<1, RAJA::loop_exec, 
    RAJA::statement::For<2, RAJA::loop_exec, 
    RAJA::statement::Lambda<0> > > > >;

  RAJA::kernel<KJI_EXECPOL>( RAJA::make_tuple(XRange, YRange, ZRange),
                             [coef0,lx,ly,lz,&hdx_2,&hdy_2,&hdz_2,
                              coefxViewConst,coefyViewConst,coefzViewConst,
                              uViewConst,etaViewConst,vpViewConst,vView,phiView] (RAJA::Index_type const i,
                                                                                  RAJA::Index_type const j,
                                                                                  RAJA::Index_type const k)
                             {
                               const float lap = (coef0*uViewConst(i+lx,j+ly,k+lz) 
                                                  +coefxViewConst(1)*(uViewConst(i+1+lx,j+ly,k+lz)+uViewConst(i-1+lx,j+ly,k+lz))
                                                  +coefyViewConst(1)*(uViewConst(i+lx,j+1+ly,k+lz)+uViewConst(i+lx,j-1+ly,k+lz))
                                                  +coefzViewConst(1)*(uViewConst(i+lx,j+ly,k+1+lz)+uViewConst(i+lx,j+ly,k-1+lz))
                                                  +coefxViewConst(2)*(uViewConst(i+2+lx,j+ly,k+lz)+uViewConst(i-2+lx,j+ly,k+lz))
                                                  +coefyViewConst(2)*(uViewConst(i+lx,j+2+ly,k+lz)+uViewConst(i+lx,j-2+ly,k+lz))
                                                  +coefzViewConst(2)*(uViewConst(i+lx,j+ly,k+2+lz)+uViewConst(i+lx,j+ly,k-2+lz))
                                                  +coefxViewConst(3)*(uViewConst(i+3+lx,j+ly,k+lz)+uViewConst(i-3+lx,j+ly,k+lz))
                                                  +coefyViewConst(3)*(uViewConst(i+lx,j+3+ly,k+lz)+uViewConst(i+lx,j-3+ly,k+lz))
                                                  +coefzViewConst(3)*(uViewConst(i+lx,j+ly,k+3+lz)+uViewConst(i+lx,j+ly,k-3+lz))
                                                  +coefxViewConst(4)*(uViewConst(i+4+lx,j+ly,k+lz)+uViewConst(i-4+lx,j+ly,k+lz))
                                                  +coefyViewConst(4)*(uViewConst(i+lx,j+4+ly,k+lz)+uViewConst(i+lx,j-4+ly,k+lz))
                                                  +coefzViewConst(4)*(uViewConst(i+lx,j+ly,k+4+lz)+uViewConst(i+lx,j+ly,k-4+lz)));
                               
                               vView(i+lx,j+ly,k+lz) =
                                 ( (2.f-etaViewConst(i+1,j+1,k+1)*etaViewConst(i+1,j+1,k+1)
                                    +2.f*etaViewConst(i+1,j+1,k+1))*uViewConst(i+lx,j+ly,k+lz)
                                   -vView(i+lx,j+ly,k+lz)
                                   +vpViewConst(i,j,k)*(lap+phiView(i,j,k)) )
                                 / (1.f+2.f*etaViewConst(i+1,j+1,k+1));
                               phiView(i,j,k) =
                                 ( phiView(i,j,k)-
                                   ( (etaViewConst(i+2,j+1,k+1)-etaViewConst(i,j+1,k+1))
                                     *(uViewConst(i+1+lx,j+ly,k+lz)-uViewConst(i-1+lx,j+ly,k+lz))*hdx_2
                                     +(etaViewConst(i+1,j+2,k+1)-etaViewConst(i+1,j,k+1))
                                     *(uViewConst(i+lx,j+1+ly,k+lz)-uViewConst(i+lx,j-1+ly,k+lz))*hdy_2
                                     +(etaViewConst(i+1,j+1,k+2)-etaViewConst(i+1,j+1,k))
                                     *(uViewConst(i+lx,j+ly,k+1+lz)-uViewConst(i+lx,j+ly,k-1+lz))*hdz_2 ) )
                                 / (1.f+etaViewConst(i+1,j+1,k+1));
                             });
#else
  for (llint i = x3; i < x4; ++i) {
       for (llint j = y3; j < y4; ++j) {
           for (llint k = z3; k < z4; ++k) {
             const float lap = (coef0*uViewConst(i+lx,j+ly,k+lz) 
                                +coefxViewConst(1)*(uViewConst(i+1+lx,j+ly,k+lz)+uViewConst(i-1+lx,j+ly,k+lz))
                                +coefyViewConst(1)*(uViewConst(i+lx,j+1+ly,k+lz)+uViewConst(i+lx,j-1+ly,k+lz))
                                +coefzViewConst(1)*(uViewConst(i+lx,j+ly,k+1+lz)+uViewConst(i+lx,j+ly,k-1+lz))
                                +coefxViewConst(2)*(uViewConst(i+2+lx,j+ly,k+lz)+uViewConst(i-2+lx,j+ly,k+lz))
                                +coefyViewConst(2)*(uViewConst(i+lx,j+2+ly,k+lz)+uViewConst(i+lx,j-2+ly,k+lz))
                                +coefzViewConst(2)*(uViewConst(i+lx,j+ly,k+2+lz)+uViewConst(i+lx,j+ly,k-2+lz))
                                +coefxViewConst(3)*(uViewConst(i+3+lx,j+ly,k+lz)+uViewConst(i-3+lx,j+ly,k+lz))
                                +coefyViewConst(3)*(uViewConst(i+lx,j+3+ly,k+lz)+uViewConst(i+lx,j-3+ly,k+lz))
                                +coefzViewConst(3)*(uViewConst(i+lx,j+ly,k+3+lz)+uViewConst(i+lx,j+ly,k-3+lz))
                                +coefxViewConst(4)*(uViewConst(i+4+lx,j+ly,k+lz)+uViewConst(i-4+lx,j+ly,k+lz))
                                +coefyViewConst(4)*(uViewConst(i+lx,j+4+ly,k+lz)+uViewConst(i+lx,j-4+ly,k+lz))
                                +coefzViewConst(4)*(uViewConst(i+lx,j+ly,k+4+lz)+uViewConst(i+lx,j+ly,k-4+lz)));
             
             vView(i+lx,j+ly,k+lz) =
               ( (2.f-etaViewConst(i+1,j+1,k+1)*etaViewConst(i+1,j+1,k+1)
                  +2.f*etaViewConst(i+1,j+1,k+1))*uViewConst(i+lx,j+ly,k+lz)
                 -vView(i+lx,j+ly,k+lz)
                 +vpViewConst(i,j,k)*(lap+phiView(i,j,k)) )
               / (1.f+2.f*etaViewConst(i+1,j+1,k+1));
             phiView(i,j,k) =
               ( phiView(i,j,k)-
                 ( (etaViewConst(i+2,j+1,k+1)-etaViewConst(i,j+1,k+1))
                   *(uViewConst(i+1+lx,j+ly,k+lz)-uViewConst(i-1+lx,j+ly,k+lz))*hdx_2
                   +(etaViewConst(i+1,j+2,k+1)-etaViewConst(i+1,j,k+1))
                   *(uViewConst(i+lx,j+1+ly,k+lz)-uViewConst(i+lx,j-1+ly,k+lz))*hdy_2
                   +(etaViewConst(i+1,j+1,k+2)-etaViewConst(i+1,j+1,k))
                   *(uViewConst(i+lx,j+ly,k+1+lz)-uViewConst(i+lx,j+ly,k-1+lz))*hdz_2 ) )
               / (1.f+etaViewConst(i+1,j+1,k+1));
             }
           }
         }
#endif
}

void target_pml_3d(llint nx, llint ny,llint nz,
                   llint x3, llint x4, llint y3,
                   llint y4, llint z3, llint z4,
                   llint lx, llint ly, llint lz,
                   float hdx_2, float hdy_2, float hdz_2,
                   const float *__restrict__ coefx,
                   const float *__restrict__ coefy,
                   const float *__restrict__ coefz,
                   const float *__restrict__ u,
                   float *__restrict__ v,
                   const float *__restrict__ vp,
                   float *__restrict__ phi,
                   const float *__restrict__ eta)
{
  const float coef0 = coefx[0] + coefy[0] + coefz[0];
  static bool haveRun = false;

  RAJA_UNUSED_VAR( nx );
  RAJA_UNUSED_VAR( coef0 );

#define KERNEL_OPTION 1
#if KERNEL_OPTION == 0
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using compute_pml_3d" << std::endl;
  }  

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxViewConst( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyViewConst( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzViewConst( coefz, 5 );

  RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > uViewConst( u, RAJA::Layout<3>(nx+2*lx, ny+2*ly, nz+2*lz) );
  RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > vView( v, RAJA::Layout<3>(nx+2*lx, ny+2*ly, nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > vpViewConst( vp, RAJA::Layout<3>(nx, ny, nz) );
  RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > phiView( phi, RAJA::Layout<3>(nx, ny, nz) );
  RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > etaViewConst( eta, RAJA::Layout<3>(nx+2, ny+2, nz+2) );  

  compute_pml_3d( RAJA::Index_type(x3), RAJA::Index_type(x4),
                  RAJA::Index_type(y3), RAJA::Index_type(y4),
                  RAJA::Index_type(z3), RAJA::Index_type(z4),
                  RAJA::Index_type(lx), RAJA::Index_type(ly), RAJA::Index_type(lz),
                  hdx_2, hdy_2, hdz_2,
                  coefxViewConst, coefyViewConst, coefzViewConst,
                  uViewConst, vpViewConst, etaViewConst, vView, phiView );
#endif
#if KERNEL_OPTION == 1
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using RAJA kernel with pointers" << std::endl;
  }

  RAJA::RangeSegment const XRange(x3, x4);
  RAJA::RangeSegment const YRange(y3, y4);
  RAJA::RangeSegment const ZRange(z3, z4);
  
  using POLICY = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::loop_exec, 
      RAJA::statement::For<1, RAJA::loop_exec, 
        RAJA::statement::For<2, RAJA::loop_exec, 
          RAJA::statement::Lambda<0>
        >
      >
    >
  >;

  RAJA::kernel<POLICY>( RAJA::make_tuple(XRange, YRange, ZRange),
    [=] ( RAJA::Index_type const i, RAJA::Index_type const j, RAJA::Index_type const k)
    {
      float const lap = LAP;

      v[IDX3_l(i,j,k)] =
         ((2.f-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)]
         +2.f*eta[IDX3_eta1(i,j,k)])*u[IDX3_l(i,j,k)]
         -v[IDX3_l(i,j,k)]
         +vp[IDX3(i,j,k)]*
         (lap+phi[IDX3(i,j,k)]))/(1.f+2.f*eta[IDX3_eta1(i,j,k)]);

      phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]-
        ((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])
             *(u[IDX3_l(i+1,j,k)]-u[IDX3_l(i-1,j,k)])*hdx_2
             +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])
             *(u[IDX3_l(i,j+1,k)]-u[IDX3_l(i,j-1,k)])*hdy_2
             +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])
           *(u[IDX3_l(i,j,k+1)]-u[IDX3_l(i,j,k-1)])*hdz_2))
             /(1.f+eta[IDX3_eta1(i,j,k)]);
    }
 );
#endif
#if KERNEL_OPTION == 2
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using original for-loops" << std::endl;
  }

    for (llint i = x3; i < x4; ++i) {
       for (llint j = y3; j < y4; ++j) {
           for (llint k = z3; k < z4; ++k) {
               float const lap = LAP;

               v[IDX3_l(i,j,k)] =
                   ((2.f-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)]
                   +2.f*eta[IDX3_eta1(i,j,k)])*u[IDX3_l(i,j,k)]
                   -v[IDX3_l(i,j,k)]
                   +vp[IDX3(i,j,k)]*
                   (lap+phi[IDX3(i,j,k)]))/(1.f+2.f*eta[IDX3_eta1(i,j,k)]);

              phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]-
                ((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])
                     *(u[IDX3_l(i+1,j,k)]-u[IDX3_l(i-1,j,k)])*hdx_2
                     +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])
                     *(u[IDX3_l(i,j+1,k)]-u[IDX3_l(i,j-1,k)])*hdy_2
                     +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])
                   *(u[IDX3_l(i,j,k+1)]-u[IDX3_l(i,j,k-1)])*hdz_2))
                     /(1.f+eta[IDX3_eta1(i,j,k)]);
          }
       }
    }
#endif
}

void target_3d(llint nx, llint ny, llint nz,
               llint x1, llint x2, llint x3,
               llint x4, llint x5, llint x6,
               llint y1, llint y2, llint y3,
               llint y4, llint y5, llint y6,
               llint z1, llint z2, llint z3,
               llint z4, llint z5, llint z6,
               llint lx, llint ly, llint lz,
               float hdx_2, float hdy_2, float hdz_2,
               const float *__restrict__ coefx,
               const float *__restrict__ coefy,
               const float *__restrict__ coefz,
               const float *__restrict__ u,
               float *__restrict__ v,
               const float *__restrict__ vp,
               float *__restrict__ phi,
               const float *__restrict__ eta)
{
    const llint xmin = 0; const llint xmax = nx;
    const llint ymin = 0; const llint ymax = ny;

    target_pml_3d(nx,ny,nz,
                  xmin,xmax,ymin,ymax,z1,z2,
                  lx,ly,lz,
                  hdx_2,hdy_2,hdz_2,
                  coefx,coefy,coefz,
                  u,v,vp,phi,eta);

    target_pml_3d(nx,ny,nz,
                  xmin,xmax,y1,y2,z3,z4,
                  lx,ly,lz,
                  hdx_2,hdy_2,hdz_2,
                  coefx,coefy,coefz,
                  u,v,vp,phi,eta);

    target_pml_3d(nx,ny,nz,
                  x1,x2,y3,y4,z3,z4,
                  lx,ly,lz,
                  hdx_2,hdy_2,hdz_2,
                  coefx,coefy,coefz,
                  u,v,vp,phi,eta);

    target_inner_3d(nx,ny,nz,
                    x3,x4,y3,y4,z3,z4,
                    lx,ly,lz,
                    coefx,coefy,coefz,
                    u,v,vp);

    target_pml_3d(nx,ny,nz,
                  x5,x6,y3,y4,z3,z4,
                  lx,ly,lz,
                  hdx_2,hdy_2,hdz_2,
                  coefx,coefy,coefz,
                  u,v,vp,phi,eta);

    target_pml_3d(nx,ny,nz,
                  xmin,xmax,y5,y6,z3,z4,
                  lx,ly,lz,
                  hdx_2,hdy_2,hdz_2,
                  coefx,coefy,coefz,
                  u,v,vp,phi,eta);

    target_pml_3d(nx,ny,nz,
                  xmin,xmax,ymin,ymax,z5,z6,
                  lx,ly,lz,
                  hdx_2,hdy_2,hdz_2,
                  coefx,coefy,coefz,
                  u,v,vp,phi,eta);
}
