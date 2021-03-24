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

#define LAPVIEW (coef0*uView(IDX3_l(i,j,k)) \
                +coefxView(1)*(uView(IDX3_l(i+1,j,k))+uView(IDX3_l(i-1,j,k))) \
                +coefyView(1)*(uView(IDX3_l(i,j+1,k))+uView(IDX3_l(i,j-1,k))) \
                +coefzView(1)*(uView(IDX3_l(i,j,k+1))+uView(IDX3_l(i,j,k-1))) \
                +coefxView(2)*(uView(IDX3_l(i+2,j,k))+uView(IDX3_l(i-2,j,k))) \
                +coefyView(2)*(uView(IDX3_l(i,j+2,k))+uView(IDX3_l(i,j-2,k))) \
                +coefzView(2)*(uView(IDX3_l(i,j,k+2))+uView(IDX3_l(i,j,k-2))) \
                +coefxView(3)*(uView(IDX3_l(i+3,j,k))+uView(IDX3_l(i-3,j,k))) \
                +coefyView(3)*(uView(IDX3_l(i,j+3,k))+uView(IDX3_l(i,j-3,k))) \
                +coefzView(3)*(uView(IDX3_l(i,j,k+3))+uView(IDX3_l(i,j,k-3))) \
                +coefxView(4)*(uView(IDX3_l(i+4,j,k))+uView(IDX3_l(i-4,j,k))) \
                +coefyView(4)*(uView(IDX3_l(i,j+4,k))+uView(IDX3_l(i,j-4,k))) \
                +coefzView(4)*(uView(IDX3_l(i,j,k+4))+uView(IDX3_l(i,j,k-4))))

#define VUPDATE v[IDX3_l(i,j,k)] = ( (2.f-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)] \
                                      +2.f*eta[IDX3_eta1(i,j,k)])*u[IDX3_l(i,j,k)]    \
                                     -v[IDX3_l(i,j,k)]                                \
                                     +vp[IDX3(i,j,k)]*(lap+phi[IDX3(i,j,k)]) )        \
                                     / (1.f+2.f*eta[IDX3_eta1(i,j,k)]);
#define VVIEWUPDATE vView(IDX3_l(i,j,k)) = ( (2.f-etaView(IDX3_eta1(i,j,k))*etaView(IDX3_eta1(i,j,k)) \
                                              +2.f*etaView(IDX3_eta1(i,j,k)))*uView(IDX3_l(i,j,k))    \
                                              -vView(IDX3_l(i,j,k))                                   \
                                              +vpView(IDX3(i,j,k))*(lap+phiView(IDX3(i,j,k))) )       \
                                              /(1.f+2.f*etaView(IDX3_eta1(i,j,k))); 

#define PHIUPDATE phi[IDX3(i,j,k)] = ( phi[IDX3(i,j,k)]                                     \
                                      - ( (eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)]) \
                                        *(u[IDX3_l(i+1,j,k)]-u[IDX3_l(i-1,j,k)])*hdx_2      \
                                        +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])  \
                                        *(u[IDX3_l(i,j+1,k)]-u[IDX3_l(i,j-1,k)])*hdy_2      \
                                        +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])  \
                                        *(u[IDX3_l(i,j,k+1)]-u[IDX3_l(i,j,k-1)])*hdz_2 ) )  \
                                        /(1.f+eta[IDX3_eta1(i,j,k)]); 
#define PHIVIEWUPDATE phiView(IDX3(i,j,k)) = ( phiView(IDX3(i,j,k))                                        \
                                            - ( (etaView(IDX3_eta1(i+1,j,k))-etaView(IDX3_eta1(i-1,j,k)))  \
                                                *(uView(IDX3_l(i+1,j,k))-uView(IDX3_l(i-1,j,k)))*hdx_2     \
                                                +(etaView(IDX3_eta1(i,j+1,k))-etaView(IDX3_eta1(i,j-1,k))) \
                                                *(uView(IDX3_l(i,j+1,k))-uView(IDX3_l(i,j-1,k)))*hdy_2     \
                                                +(etaView(IDX3_eta1(i,j,k+1))-etaView(IDX3_eta1(i,j,k-1))) \
                                                *(uView(IDX3_l(i,j,k+1))-uView(IDX3_l(i,j,k-1)))*hdz_2))   \
                                                /(1.f+etaView(IDX3_eta1(i,j,k)));


void compute_inner_3d(RAJA::Index_type const x3, const RAJA::Index_type x4, const RAJA::Index_type y3,
                      RAJA::Index_type const y4, const RAJA::Index_type z3, const RAJA::Index_type z4,
                      RAJA::Index_type const lx, const RAJA::Index_type ly, const RAJA::Index_type lz,
                      RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefxView,
                      RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefyView,
                      RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefzView,                    
                      RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ uView,
                      RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ vpView,
                      RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > &__restrict__ vView )
{
  const float coef0 = coefxView(0) + coefyView(0) + coefzView(0);

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
                                coefxView,coefyView,coefzView,
                                uView,vpView,vView] (RAJA::Index_type const i, RAJA::Index_type const j, RAJA::Index_type const k) 
                               {
                                 const float lap = (coef0*uView(i+lx,j+ly,k+lz) 
                                                    +coefxView(1)*(uView(i+1+lx,j+ly,k+lz)+uView(i-1+lx,j+ly,k+lz))
                                                    +coefyView(1)*(uView(i+lx,j+1+ly,k+lz)+uView(i+lx,j-1+ly,k+lz))
                                                    +coefzView(1)*(uView(i+lx,j+ly,k+1+lz)+uView(i+lx,j+ly,k-1+lz))
                                                    +coefxView(2)*(uView(i+2+lx,j+ly,k+lz)+uView(i-2+lx,j+ly,k+lz))
                                                    +coefyView(2)*(uView(i+lx,j+2+ly,k+lz)+uView(i+lx,j-2+ly,k+lz))
                                                    +coefzView(2)*(uView(i+lx,j+ly,k+2+lz)+uView(i+lx,j+ly,k-2+lz))
                                                    +coefxView(3)*(uView(i+3+lx,j+ly,k+lz)+uView(i-3+lx,j+ly,k+lz))
                                                    +coefyView(3)*(uView(i+lx,j+3+ly,k+lz)+uView(i+lx,j-3+ly,k+lz))
                                                    +coefzView(3)*(uView(i+lx,j+ly,k+3+lz)+uView(i+lx,j+ly,k-3+lz))
                                                    +coefxView(4)*(uView(i+4+lx,j+ly,k+lz)+uView(i-4+lx,j+ly,k+lz))
                                                    +coefyView(4)*(uView(i+lx,j+4+ly,k+lz)+uView(i+lx,j-4+ly,k+lz))
                                                    +coefzView(4)*(uView(i+lx,j+ly,k+4+lz)+uView(i+lx,j+ly,k-4+lz)));
                                 vView( i+lx,j+ly,k+lz ) =
                                   2.f*uView( i+lx,j+ly,k+lz )
                                   -vView( i+lx,j+ly,k+lz )
                                   +vpView( i,j,k )*lap;
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
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxView( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyView( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzView( coefz, 5 );
 
  RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > uView( u, RAJA::Layout<3>(nx+2*lx, ny+2*ly, nz+2*lz) );
  RAJA::View< float, RAJA::Layout<3, RAJA::Index_type, 2> > vView( v, RAJA::Layout<3>(nx+2*lx, ny+2*ly, nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<3, RAJA::Index_type, 2> > vpView( vp, RAJA::Layout<3>(nx, ny, nz) );    
    
  compute_inner_3d(RAJA::Index_type(x3), RAJA::Index_type(x4),
                   RAJA::Index_type(y3), RAJA::Index_type(y4),
                   RAJA::Index_type(z3), RAJA::Index_type(z4),
                   RAJA::Index_type(lx), RAJA::Index_type(ly), RAJA::Index_type(lz),
                   coefxView, coefyView, coefzView,
                   uView, vpView, vView );

}


void compute_pml_3d(const RAJA::Index_type ny, const RAJA::Index_type nz,
                    RAJA::Index_type const x3, const RAJA::Index_type x4, const RAJA::Index_type y3,
                    RAJA::Index_type const y4, const RAJA::Index_type z3, const RAJA::Index_type z4,
                    RAJA::Index_type const lx, const RAJA::Index_type ly, const RAJA::Index_type lz,
                    float const & hdx_2, float const & hdy_2, float const & hdz_2,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > & coefxView,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > & coefyView,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > & coefzView,                      
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > & uView,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > & vpView,
                    RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > & etaView,             
                    RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> >  vView,
                    RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> >  phiView )
{
  const float coef0 = coefxView(0) + coefyView(0) + coefzView(0);

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
      float const lap = LAPVIEW;
      VVIEWUPDATE
      PHIVIEWUPDATE     
    });

}

void compute_pml_3d_restrict(const RAJA::Index_type ny, const RAJA::Index_type nz,
                             RAJA::Index_type const x3, const RAJA::Index_type x4, const RAJA::Index_type y3,
                             RAJA::Index_type const y4, const RAJA::Index_type z3, const RAJA::Index_type z4,
                             RAJA::Index_type const lx, const RAJA::Index_type ly, const RAJA::Index_type lz,
                             float const & hdx_2, float const & hdy_2, float const & hdz_2,
                             RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefxView,
                             RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefyView,
                             RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefzView,
                             RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ uView,
                             RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ vpView,
                             RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ etaView,             
                             RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ vView,
                             RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ phiView )
{
  const float coef0 = coefxView(0) + coefyView(0) + coefzView(0);

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
      float const lap = LAPVIEW;
      VVIEWUPDATE
      PHIVIEWUPDATE     
    });
}

void compute_pml_3d_for_restrict(const RAJA::Index_type ny, const RAJA::Index_type nz,
                                 RAJA::Index_type const x3, const RAJA::Index_type x4, const RAJA::Index_type y3,
                                 RAJA::Index_type const y4, const RAJA::Index_type z3, const RAJA::Index_type z4,
                                 RAJA::Index_type const lx, const RAJA::Index_type ly, const RAJA::Index_type lz,
                                 float const & hdx_2, float const & hdy_2, float const & hdz_2,
                                 RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefxView,
                                 RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefyView,
                                 RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ coefzView,
                                 RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ uView,
                                 RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ vpView,
                                 RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ etaView,             
                                 RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ vView,
                                 RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > &__restrict__ phiView )
{
  const float coef0 = coefxView(0) + coefyView(0) + coefzView(0);

  for (llint i = x3; i < x4; ++i) {
    for (llint j = y3; j < y4; ++j) {
      for (llint k = z3; k < z4; ++k) {
        float const lap = LAPVIEW;
        VVIEWUPDATE
        PHIVIEWUPDATE   
      }
    }
  }
}


void compute_pml_3d_pointers(const RAJA::Index_type ny, const RAJA::Index_type nz,
                             RAJA::Index_type const x3, const RAJA::Index_type x4, const RAJA::Index_type y3,
                             RAJA::Index_type const y4, const RAJA::Index_type z3, const RAJA::Index_type z4,
                             RAJA::Index_type const lx, const RAJA::Index_type ly, const RAJA::Index_type lz,
                             float const & hdx_2, float const & hdy_2, float const & hdz_2,
                             const float *__restrict__ coefx,
                             const float *__restrict__ coefy,
                             const float *__restrict__ coefz,
                             const float *__restrict__ u,
                             const float *__restrict__ vp,
                             const float *__restrict__ eta,
                             float *__restrict__ v,                          
                             float *__restrict__ phi )
{
  float const coef0 = coefx[0] + coefy[0] + coefz[0];
  
  for (llint i = x3; i < x4; ++i) {
    for (llint j = y3; j < y4; ++j) {
      for (llint k = z3; k < z4; ++k) {
        float const lap = LAP;
        VUPDATE
        PHIUPDATE  
      }
    }
  }
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

#define KERNEL_OPTION 4
#if KERNEL_OPTION == 0
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using RAJA::kernel with RAJA:::View" << std::endl;
  }  

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxView( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyView( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzView( coefz, 5 );

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > uView( u, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > vView( v, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > vpView( vp, nx*ny*nz );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > phiView( phi, nx*ny*nz );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > etaView( eta, (nx+2)*(ny+2)*(nz+2) );  

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
      float const lap = LAPVIEW;
      VVIEWUPDATE
      PHIVIEWUPDATE     
    });
  
#endif
#if KERNEL_OPTION == 1
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using for loops with RAJA::View" << std::endl;
  }  

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxView( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyView( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzView( coefz, 5 );

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > uView( u, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > vView( v, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > vpView( vp, nx*ny*nz );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > phiView( phi, nx*ny*nz );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > etaView( eta, (nx+2)*(ny+2)*(nz+2) );  

  for (llint i = x3; i < x4; ++i) {
    for (llint j = y3; j < y4; ++j) {
      for (llint k = z3; k < z4; ++k) {
        float const lap = LAPVIEW;
        VVIEWUPDATE
        PHIVIEWUPDATE   
      }
    }
  }
#endif
#if KERNEL_OPTION == 2
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using RAJA::kernel inside function (with restrict) with RAJA:::View" << std::endl;
  }  

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxView( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyView( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzView( coefz, 5 );

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > uView( u, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > vView( v, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > vpView( vp, nx*ny*nz );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > phiView( phi, nx*ny*nz );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > etaView( eta, (nx+2)*(ny+2)*(nz+2) );  

  compute_pml_3d_restrict( RAJA::Index_type(ny), RAJA::Index_type(nz),
                           RAJA::Index_type(x3), RAJA::Index_type(x4),
                           RAJA::Index_type(y3), RAJA::Index_type(y4),
                           RAJA::Index_type(z3), RAJA::Index_type(z4),
                           RAJA::Index_type(lx), RAJA::Index_type(ly), RAJA::Index_type(lz),
                           hdx_2, hdy_2, hdz_2,
                           coefxView, coefyView, coefzView,
                           uView, vpView, etaView, vView, phiView ); 
#endif
#if KERNEL_OPTION == 3
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using RAJA::kernel inside function (no restrict) with RAJA:::View" << std::endl;
  }  

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxView( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyView( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzView( coefz, 5 );

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > uView( u, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > vView( v, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > vpView( vp, nx*ny*nz );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > phiView( phi, nx*ny*nz );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > etaView( eta, (nx+2)*(ny+2)*(nz+2) );  

  compute_pml_3d( RAJA::Index_type(ny), RAJA::Index_type(nz),
                  RAJA::Index_type(x3), RAJA::Index_type(x4),
                  RAJA::Index_type(y3), RAJA::Index_type(y4),
                  RAJA::Index_type(z3), RAJA::Index_type(z4),
                  RAJA::Index_type(lx), RAJA::Index_type(ly), RAJA::Index_type(lz),
                  hdx_2, hdy_2, hdz_2,
                  coefxView, coefyView, coefzView,
                  uView, vpView, etaView, vView, phiView ); 
#endif
#if KERNEL_OPTION == 4
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using for loops inside function (with restrict) with RAJA:::View" << std::endl;
  }  

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefxView( coefx, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefyView( coefy, 5 );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > coefzView( coefz, 5 );

  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > uView( u, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > vView( v, (nx+2*lx)*(ny+2*ly)*(nz+2*lz) );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > vpView( vp, nx*ny*nz );
  RAJA::View< float, RAJA::Layout<1, RAJA::Index_type, 0> > phiView( phi, nx*ny*nz );
  RAJA::View< const float, RAJA::Layout<1, RAJA::Index_type, 0> > etaView( eta, (nx+2)*(ny+2)*(nz+2) );  

  compute_pml_3d( RAJA::Index_type(ny), RAJA::Index_type(nz),
                  RAJA::Index_type(x3), RAJA::Index_type(x4),
                  RAJA::Index_type(y3), RAJA::Index_type(y4),
                  RAJA::Index_type(z3), RAJA::Index_type(z4),
                  RAJA::Index_type(lx), RAJA::Index_type(ly), RAJA::Index_type(lz),
                  hdx_2, hdy_2, hdz_2,
                  coefxView, coefyView, coefzView,
                  uView, vpView, etaView, vView, phiView ); 
#endif
#if KERNEL_OPTION == 5
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
      VUPDATE
      PHIUPDATE 
    }
 );
#endif
#if KERNEL_OPTION == 6
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using original for-loops with pointers" << std::endl;
  }

    for (llint i = x3; i < x4; ++i) {
       for (llint j = y3; j < y4; ++j) {
           for (llint k = z3; k < z4; ++k) {
               float const lap = LAP;
               VUPDATE
               PHIUPDATE        
          }
       }
    }
#endif
#if KERNEL_OPTION == 7
  if( !haveRun )
  {
    haveRun = true;
    std::cout << "Using original for-loops inside function with pointers" << std::endl;
  }

  compute_pml_3d_pointers( RAJA::Index_type(ny), RAJA::Index_type(nz),
                           RAJA::Index_type(x3), RAJA::Index_type(x4),
                           RAJA::Index_type(y3), RAJA::Index_type(y4),
                           RAJA::Index_type(z3), RAJA::Index_type(z4),
                           RAJA::Index_type(lx), RAJA::Index_type(ly), RAJA::Index_type(lz),
                           hdx_2, hdy_2, hdz_2,
                           coefx, coefy, coefz,
                           u, vp, eta, v, phi ); 

  
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
