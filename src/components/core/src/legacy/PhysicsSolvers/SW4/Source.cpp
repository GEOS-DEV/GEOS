//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "GridPointSource.h"
//#include "Source.h"
//#include "Require.h"

//#include <fenv.h>
#include <cmath>

//#include  "EW.h"
#include "Filter.h"
#include "Qspline.h"
#include "time_functions.h"
#include "Source.h"

#include "../../Common/GPException.h"

using namespace std;

//-----------------------------------------------------------------------
// Constructor, 
//
//
//    pars is only used in the 'Discrete' time function
//        when pars[1],..pars[npts] should contain the discrete function on a uniform
//        grid with spacing dt=1/freq, and pars[0] is the first time, thus the grid is
//            t_k = pars[0] + dt*k, k=0,1,..,npts-1
//    ipar should have size 1, with ipar[0] containing npts.
//
//    When the source time function is not 'Discrete', the input pars and ipars will
//    not be used.
//
//namespace SW4
//{

Source::Source( double frequency, double t0, double x0, double y0, double z0,
                double Mxx,
                double Mxy, double Mxz, double Myy, double Myz, double Mzz,
                timeDep tDep,
                const char *name, bool topodepth, int ncyc,
                double* pars,
                int npar, int* ipars, int nipar ) :
    mName( name ),
    mForces(),
    mIsMomentSource( true ),
    mFreq( frequency ),
    mT0( t0 ),
    m_myPoint( false ),
    m_zRelativeToTopography( topodepth ),
    mX0( x0 ),
    mY0( y0 ),
    mZ0( z0 ),
    mPar( nullptr ),
    mIpar( nullptr ),
    mNpar( npar ),
    mNipar( nipar ),
    mNcyc( ncyc ),
    mTimeDependence( tDep ),
    m_is_filtered( false ),
    m_zTopo( -1e38 ),
    mIgnore( false )
{
  mForces.resize( 6 );
  mForces[0] = Mxx;
  mForces[1] = Mxy;
  mForces[2] = Mxz;
  mForces[3] = Myy;
  mForces[4] = Myz;
  mForces[5] = Mzz;

  if( mNpar > 0 )
  {
    mPar = new double[mNpar];
    for( int i = 0 ; i < mNpar ; i++ )
      mPar[i] = pars[i];
  }
  else
  {
    mNpar = 2;
    mPar = new double[2];
  }

  if( mNipar > 0 )
  {
    mIpar = new int[mNipar];
    for( int i = 0 ; i < mNipar ; i++ )
      mIpar[i] = ipars[i];
  }
  else
  {
    mNipar = 1;
    mIpar = new int[1];
  }

  if( mTimeDependence == iDiscrete )
    spline_interpolation();
  else
  {
    mPar[0] = find_min_exponent();
    mPar[1] = mNcyc;
  }

// Correct source location for discrepancy between raw and smoothed topography
  correct_Z_level(); // also sets the ignore flag for sources that are above the topography

}

//-----------------------------------------------------------------------
Source::Source( double frequency, double t0, double x0, double y0, double z0,
                double Fx,
                double Fy, double Fz,
                timeDep tDep,
                const char *name, bool topodepth, int ncyc,
                double* pars,
                int npar, int* ipars, int nipar ) :
    mName( name ),
    mForces(),
    mIsMomentSource( false ),
    mFreq( frequency ),
    mT0( t0 ),
    m_myPoint( false ),
    m_zRelativeToTopography( topodepth ),
    mX0( x0 ),
    mY0( y0 ),
    mZ0( z0 ),
//  mGridPointSet(false),
    mNcyc( ncyc ),
    mTimeDependence( tDep ),
    //  m_derivative(-1),
    m_is_filtered( false ),
    //  mShearModulusFactor(correctForMu),
    m_zTopo( -1e38 ),
    mIgnore( false )
{
  mForces.resize( 3 );
  mForces[0] = Fx;
  mForces[1] = Fy;
  mForces[2] = Fz;

  mNpar = npar;
  if( mNpar > 0 )
  {
    mPar = new double[mNpar];
    for( int i = 0 ; i < mNpar ; i++ )
      mPar[i] = pars[i];
  }
  else
  {
    mNpar = 2;
    mPar = new double[2];
  }

  mNipar = nipar;
  if( mNipar > 0 )
  {
    mIpar = new int[mNipar];
    for( int i = 0 ; i < mNipar ; i++ )
      mIpar[i] = ipars[i];
  }
  else
  {
    mNipar = 1;
    mIpar = new int[1];
  }
  if( mTimeDependence == iDiscrete )
    spline_interpolation();
  else
  {
    mPar[0] = find_min_exponent();
    mPar[1] = mNcyc;
  }

// Correct source location for discrepancy between raw and smoothed topography
  correct_Z_level(); // also sets the ignore flag for sources that are above the topography

}

//-----------------------------------------------------------------------
Source::Source()
{

}

//-----------------------------------------------------------------------
Source::~Source()
{
  delete[] mPar;
  delete[] mIpar;
}

//-----------------------------------------------------------------------
double Source::getX0() const
{
  return mX0;
}

//-----------------------------------------------------------------------
double Source::getY0() const
{
  return mY0;
}

//-----------------------------------------------------------------------
double Source::getZ0() const
{
  return mZ0;
}

//-----------------------------------------------------------------------
double Source::getDepth() const
{
  return mZ0 - m_zTopo;
}

//-----------------------------------------------------------------------
double Source::getOffset() const
{
  return mT0;
}

//-----------------------------------------------------------------------
double Source::getFrequency() const
{
  return mFreq;
}

//-----------------------------------------------------------------------
void Source::setMaxFrequency( double max_freq )
{
  if( mFreq > max_freq )
    mFreq = max_freq;
}

//-----------------------------------------------------------------------
bool Source::isMomentSource() const
{
  return mIsMomentSource;
}

//-----------------------------------------------------------------------
void Source::getForces( double& fx, double& fy, double& fz ) const
                        {
  if( !mIsMomentSource )
  {
    fx = mForces[0];
    fy = mForces[1];
    fz = mForces[2];
  }
  else
    fx = fy = fz = 0;
}

//-----------------------------------------------------------------------
void Source::getMoments( double& mxx, double& mxy, double& mxz, double& myy, double& myz, double& mzz ) const
                         {
  if( mIsMomentSource )
  {
    mxx = mForces[0];
    mxy = mForces[1];
    mxz = mForces[2];
    myy = mForces[3];
    myz = mForces[4];
    mzz = mForces[5];
  }
  else
    mxx = mxy = mxz = myy = myz = mzz = 0;
}

//-----------------------------------------------------------------------
void Source::setMoments( double mxx, double mxy, double mxz, double myy, double myz, double mzz )
{
  if( mIsMomentSource )
  {

    mForces[0] = mxx;
    mForces[1] = mxy;
    mForces[2] = mxz;
    mForces[3] = myy;
    mForces[4] = myz;
    mForces[5] = mzz;
  }
  else
  {
    mForces[0] = mxx;
    mForces[1] = myy;
    mForces[2] = mzz;
  }
}

//-----------------------------------------------------------------------
double Source::getAmplitude() const
{
#define SQR(x) ((x)*(x))
  double amplitude = 0;
  if( mIsMomentSource )
  {
    double msqr = 0;
    for( int q = 0 ; q < 6 ; q++ )
      msqr += SQR( mForces[q] );
    //    amplitude = mAmp*sqrt(msqr/2.);
    msqr += SQR(mForces[1]) + SQR( mForces[2] ) + SQR( mForces[4] );
    amplitude = sqrt( 0.5 * msqr );
  }
  else
  {
    double fsqr = 0;
    for( int q = 0 ; q < 3 ; q++ )
      fsqr += SQR( mForces[q] );
    //    amplitude = mAmp*sqrt(fsqr);
    amplitude = sqrt( fsqr );
  }
  return amplitude;
#undef SQR
}

//-----------------------------------------------------------------------
ostream& operator<<( ostream& output, const Source& s )
{
  output << s.mName << ( s.isMomentSource() ? " moment" : " force" ) << " source term" << endl;
  output << "   Location (X,Y,Z) = " << s.mX0 << "," << s.mY0 << "," << s.mZ0 << " in grid no " << 0 << endl;
  output << "   Strength " << s.getAmplitude();
  output << "   t0 = " << s.mT0 << " freq = " << s.mFreq << endl;
  if( s.mIsMomentSource )
  {
    output << " Mxx Mxy Myy Mxz Myz Mzz = " << s.mForces[0] << " " << s.mForces[1] << " " << s.mForces[3] <<
           " "
           << s.mForces[2] << " " << s.mForces[4] << " " << s.mForces[5] << endl;
  }
  else
  {
    output << " Fx Fy Fz = " << s.mForces[0] << " " << s.mForces[1] << " " << s.mForces[2] << endl;
  }
  return output;
}

//-----------------------------------------------------------------------
void Source::limit_frequency( int ppw, double minvsoh )
{
  double freqlim = minvsoh / ( ppw );

  if( mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed || mTimeDependence == iDBrune ||
      mTimeDependence == iGaussian || mTimeDependence == iErf ||
      mTimeDependence == iVerySmoothBump || mTimeDependence == iSmoothWave ||
      mTimeDependence == iLiu || mTimeDependence == iC6SmoothBump )
  {
    if( mFreq > 2 * M_PI * freqlim )
      mFreq = 2 * M_PI * freqlim;
  }
  else
  {
    if( mFreq > freqlim )
      mFreq = freqlim;
  }
}

//-----------------------------------------------------------------------
double Source::compute_t0_increase( double t0_min ) const
                                    {
// Gaussian, GaussianInt=Erf, Ricker and RickerInt are all centered around mT0
  if( mTimeDependence == iGaussian || mTimeDependence == iErf )
    return t0_min + 6.0 / mFreq - mT0; // translating these by at least 6*sigma = 6/freq
  else if( mTimeDependence == iRicker || mTimeDependence == iRickerInt )
    return t0_min + 1.9 / mFreq - mT0; // 1.9 ?
  else
    return t0_min - mT0; // the rest of the time functions are zero for t<mT0
}

//-----------------------------------------------------------------------
void Source::adjust_t0( double dt0 )
{
  if( dt0 > 0 && !m_is_filtered )
    mT0 += dt0;
}

//-----------------------------------------------------------------------
double Source::dt_to_resolve( int ppw ) const
                              {
  double dt_resolved = 0;
  if( mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed || mTimeDependence == iDBrune )
  {
    const double t95 = 4.744 / mFreq;
    dt_resolved = t95 / ppw;
  }
  else
  {

  }
  return dt_resolved;
}

//-----------------------------------------------------------------------
int Source::ppw_to_resolve( double dt ) const
                            {
  int ppw = 1;
  if( mTimeDependence == iBrune || mTimeDependence == iBruneSmoothed || mTimeDependence == iDBrune )
  {
    const double t95 = 4.744 / mFreq;
    ppw = static_cast<int>( t95 / dt );
  }
  else
  {

  }
  return ppw;
}

//-----------------------------------------------------------------------
void Source::setFrequency( double freq )
{
  mFreq = freq;
}

//-----------------------------------------------------------------------
void Source::correct_Z_level()
{
// this routine 
// 1. calculates the z-coordinate of the topography right above the source and saves it in m_zTopo
// 2. if m_relativeToTopography == true, it adds m_zTopo to mZ0
// 3. checks if the source is inside the computational domain. If not, set mIgnore=true

//   right now, only Cartesian grid supported
//  if( !a_ew->topographyExists() ) // this is the easy case w/o topography
  {
    m_zTopo = 0.;
// make sure the station is below or on the topography (z is positive downwards)
    if( mZ0 < m_zTopo - 1e-9 )
    {
      mIgnore = true;
      cout << "Ignoring Source at X= " << mX0 << " Y= " << mY0 << "  Z= " << mZ0 <<
           ", because it is above the topography z= "
           << m_zTopo << endl;
    }
    return; // done with the flat case
  }
}

//-----------------------------------------------------------------------
void Source::getsourcewgh( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
                           {
  // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
  double p5 = ai * ai * ai * ai * ai * ( 5.0 / 3 - 7.0 / 24 * ai - 17 / 12.0 * ai * ai + 1.125 * ai * ai * ai - 0.25 * ai * ai * ai * ai );
  wgh[0] = 1.0 / 24 * ( 2 * ai - ai * ai - 2 * ai * ai * ai - 19 * ai * ai * ai * ai ) + p5;
  wgh[1] = 1.0 / 6 * ( -4 * ai + 4 * ai * ai + ai * ai * ai ) + 4 * ai * ai * ai * ai - 5 * p5;
  wgh[2] = 1 - 1.25 * ai * ai - 97.0 / 12 * ai * ai * ai * ai + 10 * p5;
  wgh[3] = 1.0 / 6 * ( 4 * ai + 4 * ai * ai - ai * ai * ai + 49 * ai * ai * ai * ai ) - 10 * p5;
  wgh[4] = 1.0 / 24 * ( -2 * ai - ai * ai + 2 * ai * ai * ai ) - 4.125 * ai * ai * ai * ai + 5 * p5;
  wgh[5] = 5.0 / 6 * ai * ai * ai * ai - p5;

  // Derivatives of wgh wrt. ai:
  p5 = 5 * ai * ai * ai * ai * ( 5.0 / 3 - 7.0 / 24 * ai - 17 / 12.0 * ai * ai + 1.125 * ai * ai * ai - 0.25 * ai * ai * ai * ai ) +
      ai * ai * ai * ai * ai * ( -7.0 / 24 - 17 / 6.0 * ai + 3 * 1.125 * ai * ai - ai * ai * ai );
  dwghda[0] = 1.0 / 24 * ( 2 - 2 * ai - 6 * ai * ai - 19 * 4 * ai * ai * ai ) + p5;
  dwghda[1] = 1.0 / 6 * ( -4 + 8 * ai + 3 * ai * ai ) + 16 * ai * ai * ai - 5 * p5;
  dwghda[2] = -2.5 * ai - 97.0 / 3 * ai * ai * ai + 10 * p5;
  dwghda[3] = 1.0 / 6 * ( 4 + 8 * ai - 3 * ai * ai + 49 * 4 * ai * ai * ai ) - 10 * p5;
  dwghda[4] = 1.0 / 24 * ( -2 - 2 * ai + 6 * ai * ai ) - 4.125 * 4 * ai * ai * ai + 5 * p5;
  dwghda[5] = 20.0 / 6 * ai * ai * ai - p5;

  // Second derivatives of wgh wrt. ai:
  p5 = ai * ai * ai * ( 100.0 / 3 - 8.75 * ai - 59.5 * ai * ai + 63 * ai * ai * ai - 18 * ai * ai * ai * ai );

  ddwghda[0] = -1.0 / 12 - 0.5 * ai - 9.5 * ai * ai + p5;
  ddwghda[1] = 4.0 / 3 + ai + 48 * ai * ai - 5 * p5;
  ddwghda[2] = -2.5 - 97 * ai * ai + 10 * p5;
  ddwghda[3] = 4.0 / 3 - ai + 98 * ai * ai - 10 * p5;
  ddwghda[4] = -1.0 / 12 + 0.5 * ai - 49.5 * ai * ai + 5 * p5;
  ddwghda[5] = 10 * ai * ai - p5;

}

//-----------------------------------------------------------------------
void Source::getsourcedwgh( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
                            {
  // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position
  double p5 = ai * ai * ai * ai * ( -25.0 / 12 - 0.75 * ai + 59.0 / 12 * ai * ai - 4 * ai * ai * ai + ai * ai * ai * ai );
  wgh[0] = 1.0 / 12 * ( -1 + ai + 3 * ai * ai + 8 * ai * ai * ai ) + p5;
  wgh[1] = 2.0 / 3 * ( 1 - 2 * ai ) - 0.5 * ai * ai - 3.5 * ai * ai * ai - 5 * p5;
  wgh[2] = 2.5 * ai + 22.0 / 3 * ai * ai * ai + 10 * p5;
  wgh[3] = 2.0 / 3 * ( -1 - 2 * ai ) + 0.5 * ai * ai - 23.0 / 3 * ai * ai * ai - 10 * p5;
  wgh[4] = ( 1 + ai ) / 12 - 0.25 * ai * ai + 4 * ai * ai * ai + 5 * p5;
  wgh[5] = -5.0 / 6 * ai * ai * ai - p5;

  // Derivatives of wgh wrt. ai:
  p5 = 4 * ai * ai * ai * ( -25.0 / 12 - 0.75 * ai + 59.0 / 12 * ai * ai - 4 * ai * ai * ai + ai * ai * ai * ai ) +
      ai * ai * ai * ai * ( -0.75 + 59.0 / 6 * ai - 12 * ai * ai + 4 * ai * ai * ai );
  dwghda[0] = 1.0 / 12 * ( 1 + 6 * ai + 24 * ai * ai ) + p5;
  dwghda[1] = 2.0 / 3 * ( -2 ) - ai - 3 * 3.5 * ai * ai - 5 * p5;
  dwghda[2] = 2.5 + 22.0 * ai * ai + 10 * p5;
  dwghda[3] = 2.0 / 3 * ( -2 ) + ai - 23.0 * ai * ai - 10 * p5;
  dwghda[4] = 1.0 / 12 - 0.5 * ai + 12 * ai * ai + 5 * p5;
  dwghda[5] = -5.0 / 2 * ai * ai - p5;

  // Second derivatives of wgh wrt. ai:
  p5 = ai * ai * ( -25 - 15 * ai + 147.5 * ai * ai - 168 * ai * ai * ai + 56 * ai * ai * ai * ai );

  ddwghda[0] = 0.5 + 4 * ai + p5;
  ddwghda[1] = -1 - 21 * ai - 5 * p5;
  ddwghda[2] = 44 * ai + 10 * p5;
  ddwghda[3] = 1 - 46 * ai - 10 * p5;
  ddwghda[4] = -0.5 + 24 * ai + 5 * p5;
  ddwghda[5] = -5 * ai - p5;

}

//-----------------------------------------------------------------------
void Source::getsourcewghlow( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
                              {
  // Lower component stencil, to use at lower boundaries
  // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position

  wgh[0] = ( 2 * ai - ai * ai - 2 * ai * ai * ai + ai * ai * ai * ai ) / 24;
  wgh[1] = ( -4 * ai + 4 * ai * ai + ai * ai * ai - ai * ai * ai * ai ) / 6;
  wgh[2] = 1 - 1.25 * ai * ai + 0.25 * ai * ai * ai * ai;
  wgh[3] = ( 4 * ai + 4 * ai * ai - ai * ai * ai - ai * ai * ai * ai ) / 6;
  wgh[4] = ( -2 * ai - ai * ai + 2 * ai * ai * ai + ai * ai * ai * ai ) / 24;
  wgh[5] = 0;

  // Derivatives of wgh wrt. ai:
  dwghda[0] = ( 1 - ai - 3 * ai * ai + 2 * ai * ai * ai ) / 12;
  dwghda[1] = ( -2 + 4 * ai + 1.5 * ai * ai - 2 * ai * ai * ai ) / 3;
  dwghda[2] = -2.5 * ai + ai * ai * ai;
  dwghda[3] = ( 2 + 4 * ai - 1.5 * ai * ai - 2 * ai * ai * ai ) / 3;
  dwghda[4] = ( -1 - ai + 3 * ai * ai + 2 * ai * ai * ai ) / 12;
  dwghda[5] = 0;

  // Second derivatives of wgh wrt. ai:
  ddwghda[0] = -1.0 / 12 - 0.5 * ai + 0.5 * ai * ai;
  ddwghda[1] = 4.0 / 3 + ai - 2 * ai * ai;
  ddwghda[2] = -2.5 + 3 * ai * ai;
  ddwghda[3] = 4.0 / 3 - ai - 2 * ai * ai;
  ddwghda[4] = -1.0 / 12 + 0.5 * ai + 0.5 * ai * ai;
  ddwghda[5] = 0;
}

//-----------------------------------------------------------------------
void Source::getsourcedwghlow( double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const
                               {
  // Lower component stencil, to use at lower boundaries, dirac derivative weights.
  // Moments k=0,1,2,3,4 exact, two cont. derivatives wrt. position

  // same as derivatives of dirac weights.
  wgh[0] = ( -1 + ai + 3 * ai * ai - 2 * ai * ai * ai ) / 12;
  wgh[1] = ( 2 - 4 * ai - 1.5 * ai * ai + 2 * ai * ai * ai ) / 3;
  wgh[2] = 2.5 * ai - ai * ai * ai;
  wgh[3] = ( -2 - 4 * ai + 1.5 * ai * ai + 2 * ai * ai * ai ) / 3;
  wgh[4] = ( 1 + ai - 3 * ai * ai - 2 * ai * ai * ai ) / 12;
  wgh[5] = 0;

  // Derivatives of wgh wrt. ai:
  dwghda[0] = 1.0 / 12 + 0.5 * ai - 0.5 * ai * ai;
  dwghda[1] = -4.0 / 3 - ai + 2 * ai * ai;
  dwghda[2] = 2.5 - 3 * ai * ai;
  dwghda[3] = -4.0 / 3 + ai + 2 * ai * ai;
  dwghda[4] = 1.0 / 12 - 0.5 * ai - 0.5 * ai * ai;
  dwghda[5] = 0;

  // Second derivatives of wgh wrt. ai:
  ddwghda[0] = 0.5 - ai;
  ddwghda[1] = -1 + 4 * ai;
  ddwghda[2] = -6 * ai;
  ddwghda[3] = 1 + 4 * ai;
  ddwghda[4] = -0.5 - ai;
  ddwghda[5] = 0;
}

//-----------------------------------------------------------------------
void Source::getmetwgh( double ai, double wgh[8], double dwgh[8],
                        double ddwgh[8],
                        double dddwgh[8] ) const
                        {
  double pol = ai * ai * ai * ai * ai * ai * ai * ( -251 + 135 * ai + 25 * ai * ai -
      33 * ai * ai * ai + 6 * ai * ai * ai * ai ) / 720;

  wgh[0] = -1.0 / 60 * ai + 1.0 / 180 * ai * ai + 1.0 / 48 * ai * ai * ai + 23.0 / 144 * ai * ai * ai * ai
      - ( 17.0 * ai + 223.0 ) * ai * ai * ai * ai * ai / 720 - pol;
  wgh[1] = 3.0 / 20 * ai - 3.0 / 40 * ai * ai - 1.0 / 6 * ai * ai * ai - 13.0 / 12 * ai * ai * ai * ai +
      97.0 / 45 * ai * ai * ai * ai * ai + 1.0 / 6 * ai * ai * ai * ai * ai * ai + 7 * pol;
  wgh[2] = -0.75 * ai + 0.75 * ai * ai + ( 13.0 + 155 * ai ) * ai * ai * ai / 48 - 103.0 / 16 * ai * ai * ai * ai * ai
      - 121.0 / 240 * ai * ai * ai * ai * ai * ai - 21 * pol;
  wgh[3] = 1 - 49.0 / 36 * ai * ai - 49.0 / 9 * ai * ai * ai * ai + 385.0 / 36 * ai * ai * ai * ai * ai +
      61.0 / 72 * ai * ai * ai * ai * ai * ai + 35 * pol;
  wgh[4] = 0.75 * ai + 0.75 * ai * ai - 13.0 / 48 * ai * ai * ai + 89.0 / 16 * ai * ai * ai * ai -
      1537.0 / 144 * ai * ai * ai * ai * ai - 41.0 / 48 * ai * ai * ai * ai * ai * ai - 35 * pol;
  wgh[5] = -3.0 / 20 * ai - 3.0 / 40 * ai * ai + 1.0 / 6 * ai * ai * ai - 41.0 / 12 * ai * ai * ai * ai
      + 6.4 * ai * ai * ai * ai * ai + 31.0 / 60 * ai * ai * ai * ai * ai * ai + 21 * pol;
  wgh[6] = 1.0 / 60 * ai + 1.0 / 180 * ai * ai - 1.0 / 48 * ai * ai * ai + 167.0 / 144 * ai * ai * ai * ai -
      1537.0 / 720 * ai * ai * ai * ai * ai - 25.0 / 144 * ai * ai * ai * ai * ai * ai - 7 * pol;
  wgh[7] = -1.0 / 6 * ai * ai * ai * ai + 11.0 / 36 * ai * ai * ai * ai * ai + 1.0 / 40 * ai * ai * ai * ai * ai * ai + pol;

  // Derivative wrt. ai
  pol = ai * ai * ai * ai * ai * ai * ( -1757.0 / 720 + 1.5 * ai + 0.31250 * ai * ai - ( 1.375 * ai * ai * ai - 0.275 * ai * ai * ai * ai ) / 3 );
  dwgh[0] = -1.0 / 60 + 1.0 / 90 * ai + ai * ai / 16 + 23.0 / 36 * ai * ai * ai - 223.0 / 144 * ai * ai * ai * ai -
      17.0 / 120 * ai * ai * ai * ai * ai - pol;
  dwgh[1] = 3.0 / 20 - 3.0 / 20 * ai - 0.5 * ai * ai - 13.0 / 3 * ai * ai * ai + 97.0 / 9 * ai * ai * ai * ai +
      ai * ai * ai * ai * ai + 7 * pol;
  dwgh[2] = -0.75 + 1.5 * ai + 13.0 / 16 * ai * ai + 155.0 * ai * ai * ai / 12 - 103.0 * 5.0 / 16 * ai * ai * ai * ai
      - 121.0 / 40 * ai * ai * ai * ai * ai - 21 * pol;
  dwgh[3] = -49.0 / 18 * ai - 4 * 49.0 / 9.0 * ai * ai * ai + 385.0 * 5.0 / 36 * ai * ai * ai * ai +
      61.0 / 12 * ai * ai * ai * ai * ai + 35 * pol;
  dwgh[4] = 0.75 + 1.5 * ai - 13.0 / 16 * ai * ai + 89.0 / 4 * ai * ai * ai - 1537.0 * 5 / 144.0 * ai * ai * ai * ai -
      41.0 / 8 * ai * ai * ai * ai * ai - 35 * pol;
  dwgh[5] = -3.0 / 20 - 3.0 / 20 * ai + 0.5 * ai * ai - 41.0 / 3 * ai * ai * ai + 32 * ai * ai * ai * ai +
      3.1 * ai * ai * ai * ai * ai + 21 * pol;
  dwgh[6] = 1.0 / 60 + 1.0 / 90 * ai - 1.0 / 16 * ai * ai + 167.0 / 36 * ai * ai * ai - 1537.0 / 144 * ai * ai * ai * ai -
      25.0 / 24 * ai * ai * ai * ai * ai - 7 * pol;
  dwgh[7] = -2.0 / 3 * ai * ai * ai + 55.0 / 36 * ai * ai * ai * ai + 3.0 / 20 * ai * ai * ai * ai * ai + pol;

  // Second derivative wrt. ai
  pol = ai * ai * ai * ai * ai * ( -1757.0 / 120 + 10.5 * ai + 2.5 * ai * ai - 4.125 * ai * ai * ai + 11.0 / 12 * ai * ai * ai * ai );
  ddwgh[0] = 1.0 / 90 + 0.125 * ai + 23.0 / 12 * ai * ai - 223.0 / 36 * ai * ai * ai - 17.0 / 24 * ai * ai * ai * ai - pol;
  ddwgh[1] = -3.0 / 20 - ai - 13.0 * ai * ai + 4 * 97.0 / 9.0 * ai * ai * ai + 5 * ai * ai * ai * ai + 7 * pol;
  ddwgh[2] = 1.5 + 13.0 / 8 * ai + 155.0 / 4 * ai * ai - 103.0 * 5.0 / 4 * ai * ai * ai - 121.0 / 8 * ai * ai * ai * ai - 21 * pol;
  ddwgh[3] = -49.0 / 18 - 4 * 49.0 / 3.0 * ai * ai + 385.0 * 5.0 / 9.0 * ai * ai * ai + 5 * 61.0 / 12 * ai * ai * ai * ai + 35 * pol;
  ddwgh[4] = 1.5 - 13.0 / 8 * ai + 89.0 * 3.0 / 4 * ai * ai - 1537.0 * 5.0 / 36 * ai * ai * ai - 205.0 / 8 * ai * ai * ai * ai - 35 * pol;
  ddwgh[5] = -3.0 / 20 + ai - 41.0 * ai * ai + 128 * ai * ai * ai + 15.5 * ai * ai * ai * ai + 21 * pol;
  ddwgh[6] = 1.0 / 90 - 0.125 * ai + 167.0 / 12 * ai * ai - 1537.0 / 36 * ai * ai * ai - 125.0 / 24 * ai * ai * ai * ai - 7 * pol;
  ddwgh[7] = -2 * ai * ai + 220.0 / 36 * ai * ai * ai + 0.75 * ai * ai * ai * ai + pol;

  // Third derivative wrt. ai
  pol = ai * ai * ai * ai * ( -1757.0 / 24 + 63 * ai + 17.5 * ai * ai - 33 * ai * ai * ai + 8.25 * ai * ai * ai * ai );
  dddwgh[0] = 0.125 + 23.0 / 6 * ai - 223.0 / 12 * ai * ai - 17.0 / 6 * ai * ai * ai - pol;
  dddwgh[1] = -1 - 26.0 * ai + 4 * 97.0 / 3 * ai * ai + 20 * ai * ai * ai + 7 * pol;
  dddwgh[2] = 1.625 + 77.5 * ai - 386.25 * ai * ai - 60.5 * ai * ai * ai - 21 * pol;
  dddwgh[3] = -392.0 / 3 * ai + 1925.0 / 3 * ai * ai + 305.0 / 3 * ai * ai * ai + 35 * pol;
  dddwgh[4] = -1.625 + 133.5 * ai - 7685.0 / 12 * ai * ai - 102.5 * ai * ai * ai - 35 * pol;
  dddwgh[5] = 1 - 82.0 * ai + 384.0 * ai * ai + 62.0 * ai * ai * ai + 21 * pol;
  dddwgh[6] = -0.125 + 167.0 / 6 * ai - 1537.0 / 12 * ai * ai - 125.0 / 6 * ai * ai * ai - 7 * pol;
  dddwgh[7] = -4 * ai + 220.0 / 12 * ai * ai + 3 * ai * ai * ai + pol;
}

//-----------------------------------------------------------------------
void Source::getmetdwgh( double ai, double wgh[8] ) const
                         {
  double pol = ai * ai * ai * ai * ai * ai * ( -827 + 420 * ai + 165 * ai * ai - 180 * ai * ai * ai
      + 36 * ai * ai * ai * ai ) / 720;

  wgh[0] = -1.0 / 60 + 1.0 / 90 * ai + 1.0 / 16 * ai * ai + 5.0 / 36 * ai * ai * ai -
      55.0 / 144 * ai * ai * ai * ai - 7.0 / 20 * ai * ai * ai * ai * ai - pol;
  wgh[1] = 3.0 / 20 * ( 1 - ai ) - 0.5 * ai * ai - 5.0 / 6 * ai * ai * ai + 47.0 / 18 * ai * ai * ai * ai
      + 59.0 / 24 * ai * ai * ai * ai * ai + 7 * pol;
  wgh[2] = -0.75 + 1.5 * ai + 13.0 / 16 * ai * ai + 29.0 / 12 * ai * ai * ai - 123.0 / 16 * ai * ai * ai * ai -
      7.4 * ai * ai * ai * ai * ai - 21 * pol;
  wgh[3] = ( -49.0 * ai - 77.0 * ai * ai * ai ) / 18 + 455.0 / 36 * ai * ai * ai * ai + 99.0 / 8 * ai * ai * ai * ai * ai
      + 35 * pol;
  wgh[4] = 0.75 + 1.5 * ai - 13.0 / 16 * ai * ai + 4.75 * ai * ai * ai - 1805.0 / 144 * ai * ai * ai * ai -
      149.0 / 12 * ai * ai * ai * ai * ai - 35 * pol;
  wgh[5] = -3.0 / 20 * ( 1 + ai ) + 0.5 * ai * ai - 19.0 / 6 * ai * ai * ai + 7.5 * ai * ai * ai * ai +
      299.0 / 40 * ai * ai * ai * ai * ai + 21 * pol;
  wgh[6] = 1.0 / 60 + 1.0 / 90 * ai - 1.0 / 16 * ai * ai + 41.0 / 36 * ai * ai * ai - 361.0 / 144 * ai * ai * ai * ai -
      2.5 * ai * ai * ai * ai * ai - 7 * pol;
  wgh[7] = -1.0 / 6 * ai * ai * ai + 13.0 / 36 * ai * ai * ai * ai + 43.0 / 120 * ai * ai * ai * ai * ai + pol;
}
//-----------------------------------------------------------------------
void Source::getmetwgh7( double ai, double wgh[7] ) const
                         {
  wgh[0] = ai * ai / 180.0 - ai * ai * ai * ai / 144.0 + ai * ai * ai * ai * ai * ai / 720.0 - ai / 60.0 +
      ai * ai * ai / 48.0 - ai * ai * ai * ai * ai / 240.0;
  wgh[1] = -3.0 / 40.0 * ai * ai + ai * ai * ai * ai / 12.0 - ai * ai * ai * ai * ai * ai / 120.0 +
      3.0 / 20.0 * ai - ai * ai * ai / 6.0 + ai * ai * ai * ai * ai / 60.0;
  wgh[2] = 3.0 / 4.0 * ai * ai - 13.0 / 48.0 * ai * ai * ai * ai + ai * ai * ai * ai * ai * ai / 48.0 -
      3.0 / 4.0 * ai + 13.0 / 48.0 * ai * ai * ai - ai * ai * ai * ai * ai / 48.0;
  wgh[3] = 1.0 - 49.0 / 36.0 * ai * ai + 7.0 / 18.0 * ai * ai * ai * ai - ai * ai * ai * ai * ai * ai / 36.0;
  wgh[4] = 3.0 / 4.0 * ai * ai - 13.0 / 48.0 * ai * ai * ai * ai + 3.0 / 4.0 * ai - 13.0 / 48.0 * ai * ai * ai +
      ai * ai * ai * ai * ai / 48.0 + ai * ai * ai * ai * ai * ai / 48.0;
  wgh[5] = -3.0 / 40.0 * ai * ai + ai * ai * ai * ai / 12.0 - 3.0 / 20.0 * ai + ai * ai * ai / 6.0 -
      ai * ai * ai * ai * ai / 60.0 - ai * ai * ai * ai * ai * ai / 120.0;
  wgh[6] = ai * ai / 180.0 - ai * ai * ai * ai / 144.0 + ai * ai * ai * ai * ai * ai / 720.0 + ai / 60.0 +
      ai * ai * ai * ai * ai / 240.0 - ai * ai * ai / 48.0;

}

//-----------------------------------------------------------------------
void Source::getmetdwgh7( double ai, double wgh[7] ) const
                          {
  wgh[0] = -1.0 / 60.0 + ai * ai / 16.0 - ai * ai * ai * ai / 48.0 + ai / 90.0 -
      ai * ai * ai / 36.0 + ai * ai * ai * ai * ai / 120.0;
  wgh[1] = 3.0 / 20.0 - ai * ai / 2.0 + ai * ai * ai * ai / 12.0 - 3.0 / 20.0 * ai +
      ai * ai * ai / 3.0 - ai * ai * ai * ai * ai / 20.0;
  wgh[2] = -3.0 / 4.0 + 13.0 / 16.0 * ai * ai - 5.0 / 48.0 * ai * ai * ai * ai +
      3.0 / 2.0 * ai - 13.0 / 12.0 * ai * ai * ai + ai * ai * ai * ai * ai / 8.0;
  wgh[3] = -49.0 / 18.0 * ai + 14.0 / 9.0 * ai * ai * ai - ai * ai * ai * ai * ai / 6.0;
  wgh[4] = 3.0 / 4.0 - 13.0 / 16.0 * ai * ai + 3.0 / 2.0 * ai - 13.0 / 12.0 * ai * ai * ai +
      5.0 / 48.0 * ai * ai * ai * ai + ai * ai * ai * ai * ai / 8.0;
  wgh[5] = -3.0 / 20.0 + ai * ai / 2.0 - ai * ai * ai * ai / 12.0 - 3.0 / 20.0 * ai +
      ai * ai * ai / 3.0 - ai * ai * ai * ai * ai / 20.0;
  wgh[6] = 1.0 / 60.0 - ai * ai / 16.0 + ai * ai * ai * ai / 48.0 + ai / 90.0 -
      ai * ai * ai / 36.0 + ai * ai * ai * ai * ai / 120.0;
}

//-----------------------------------------------------------------------
void Source::set_grid_point_sources4( double h, double xmin, double ymin,
                                      double zmin,
                                      int Ni, int Nj,
                                      int Nz,
                                      vector<GridPointSource*>& point_sources,
                                      int interior[6] ) const
                                      {
// note that this routine is called from all processors, for each input source 
//   int i,j,k,g;
  double q, r, s;
  //   double h = a_EW->mGridSize[g];
  bool canBeInverted, curvilinear;
  double normwgh[4] =
  { 17.0 / 48.0, 59.0 / 48.0, 43.0 / 48.0, 49.0 / 48.0 };
  int g = 0;

  // currently, only Cartesian grids
  {
// Cartesian case
    q = ( mX0 - xmin ) / h + 1;
    r = ( mY0 - ymin ) / h + 1;
    s = ( mZ0 - zmin ) / h + 1;
    canBeInverted = true;
    curvilinear = false;
  }

  //   int Ni = a_EW->m_global_nx[g];
  //   int Nj = a_EW->m_global_ny[g];
  //   int Nz = a_EW->m_global_nz[g];

  int ic, jc, kc;
  bool ccbndry; //upperbndry, lowerbndry
  double ai, bi, ci;

// Delta distribution
  double wghi[6], wghj[6], wghk[6], wghix[6], wghjy[6], wghkz[6];
  double wghixx[6], wghjyy[6], wghkzz[6];
// Delta' distribution
  double dwghi[6], dwghj[6], dwghk[6], dwghix[6], dwghjy[6], dwghkz[6];
  double dwghixx[6], dwghjyy[6], dwghkzz[6];

  if( canBeInverted )
  {
    // Compute source location and weights in source discretization
    ic = static_cast<int>( floor( q ) );
    jc = static_cast<int>( floor( r ) );
    kc = static_cast<int>( floor( s ) );

// Bias stencil away from boundary, no source at ghost/padding points
    if( ic <= 2 )
      ic = 3;
    if( ic >= Ni - 2 )
      ic = Ni - 3;
    if( jc <= 2 )
      jc = 3;
    if( jc >= Nj - 2 )
      jc = Nj - 3;

// Six point stencil, with points, kc-2,..kc+3, Interior in domain
// if kc-2>=1, kc+3 <= Nz --> kc >= 3 and kc <= Nz-3
// Can evaluate with two ghost points if kc-2>=-1 and kc+3 <= Nz+2
//  --->  kc >=1 and kc <= Nz-1
//
    if( kc >= Nz )
      kc = Nz - 1;
    if( kc < 1 )
      kc = 1;

// upper(surface) and lower boundaries , when the six point stencil kc-2,..kc+3
// make use of the first (k=1) or the last (k=Nz) interior point.
    //upperbndry = (kc == 1    || kc == 2    || kc == 3  );
    //lowerbndry = (kc == Nz-1 || kc == Nz-2 || kc == Nz-3 );

// ccbndry=true if at the interface between the curvilinear grid and the cartesian grid. 
// Defined as the six point stencil uses values from both grids.
//      ccbndry = a_EW->topographyExists() &&  ( (upperbndry && g == a_EW->mNumberOfGrids-2) ||
//					       (lowerbndry && g == a_EW->mNumberOfGrids-1)    );
    ccbndry = false;
// If not at the interface between curvilinear and Cartesian grids, bias stencil away 
// from the boundary.
    if( !ccbndry )
    {
      if( kc <= 2 )
        kc = 3;
      if( kc >= Nz - 2 )
        kc = Nz - 3;
    }
    ai = q - ic, bi = r - jc, ci = s - kc;

    // Delta distribution
    getsourcewgh( ai, wghi, wghix, wghixx );
    getsourcewgh( bi, wghj, wghjy, wghjyy );
    getsourcewgh( ci, wghk, wghkz, wghkzz );

// Delta' distribution
    getsourcedwgh( ai, dwghi, dwghix, dwghixx );
    getsourcedwgh( bi, dwghj, dwghjy, dwghjyy );
    getsourcedwgh( ci, dwghk, dwghkz, dwghkzz );

    // Special boundary stencil at free surface
    if( !ccbndry && ( kc == 3 && ci <= 0 ) )
    {
      getsourcewghlow( ci, wghk, wghkz, wghkzz );
      getsourcedwghlow( ci, dwghk, dwghkz, dwghkzz );
    }

// Boundary correction, at upper boundary, but only if SBP operators are used there
//
//CHANGE      if( (g == a_EW->mNumberOfGrids-1) && a_EW->is_onesided(g,4)  )
// can do without the if-statement, since there is only one grid, and we
// assume that one sided operators are always used on the top boundary.
    {
      for( int k = 0 ; k <= 5 ; k++ )
      {
        if( ( 1 <= k + kc - 2 ) && ( k + kc - 2 <= 4 ) )
        {
          wghk[k] /= normwgh[k + kc - 3];
          dwghk[k] /= normwgh[k + kc - 3];
          wghkz[k] /= normwgh[k + kc - 3];
          dwghkz[k] /= normwgh[k + kc - 3];
          wghkzz[k] /= normwgh[k + kc - 3];
          dwghkzz[k] /= normwgh[k + kc - 3];
        }
      }
    }
  }
  if( !mIsMomentSource && canBeInverted )
  {
    for( int k = kc - 2 ; k <= kc + 3 ; k++ )
      for( int j = jc - 2 ; j <= jc + 3 ; j++ )
        for( int i = ic - 2 ; i <= ic + 3 ; i++ )
        {
          double wF = wghi[i - ic + 2] * wghj[j - jc + 2] * wghk[k - kc + 2];
          if( ( wF != 0 ) && ( mForces[0] != 0 || mForces[1] != 0 || mForces[2] != 0 ) && is_inside( i, j, k,
                                                                                                     interior ) )
          //		   && a_EW->interior_point_in_proc(i,j,g) ) // checks if (i,j) belongs to this processor
          {
            //		  if( curvilinear )
            //		     wF /= a_EW->mJ(i,j,k);
            //		  else
            wF /= h * h * h;

            if( 1 <= k && k <= Nz )
            {
              GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0,
                                                                i,
                                                                j, k, g,
                                                                wF * mForces[0],
                                                                wF * mForces[1], wF * mForces[2],
                                                                mTimeDependence,
                                                                mNcyc,
                                                                mPar,
                                                                mNpar, mIpar, mNipar );
              point_sources.push_back( sourcePtr );
            }
          }
        }
  }
  // Moment source.
  else if( mIsMomentSource )
  {
    double qX0[3], rX0[3], sX0[3];
    if( !curvilinear )
    {
      // Cartesian case, constant metric
      qX0[0] = 1 / h;
      qX0[1] = 0;
      qX0[2] = 0;
      rX0[0] = 0;
      rX0[1] = 1 / h;
      rX0[2] = 0;
      sX0[0] = 0;
      sX0[1] = 0;
      sX0[2] = 1 / h;
    }
    if( canBeInverted )
    {
//         cout << "Source: h =  " << h << endl;
      for( int k = kc - 2 ; k <= kc + 3 ; k++ )
        for( int j = jc - 2 ; j <= jc + 3 ; j++ )
          for( int i = ic - 2 ; i <= ic + 3 ; i++ )
          {
            double wFx = 0, wFy = 0, wFz = 0;
            //CHANGE		  if( a_EW->interior_point_in_proc(i,j,g) )
            if( is_inside( i, j, k, interior ) )
            {
              //                     cout << " src at " << i << " " << j << " " << k << endl;
              wFx += qX0[0] * dwghi[i - ic + 2] * wghj[j - jc + 2] * wghk[k - kc + 2];
              //		  wFy += qX0[1]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2];
              //		  wFz += qX0[2]*dwghi[i-ic+2]* wghj[j-jc+2]* wghk[k-kc+2];

              //		  wFx +=  wghi[i-ic+2]*rX0[0]*dwghj[j-jc+2]* wghk[k-kc+2];
              wFy += wghi[i - ic + 2] * rX0[1] * dwghj[j - jc + 2] * wghk[k - kc + 2];
              //		  wFz +=  wghi[i-ic+2]*rX0[2]*dwghj[j-jc+2]* wghk[k-kc+2];

              //		     wFx +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[0]*dwghk[k-kc+2];
              //		     wFx +=  sX0[0];
              //		     wFy +=  wghi[i-ic+2]* wghj[j-jc+2]*sX0[1]*dwghk[k-kc+2];
              wFz += wghi[i - ic + 2] * wghj[j - jc + 2] * sX0[2] * dwghk[k - kc + 2];

              double jaci = 1.0 / ( h * h * h );
              //		     if( curvilinear )
              //			jaci = 1/a_EW->mJ(i,j,k);
              //		     else
              //			jaci =
              double fx = -( mForces[0] * wFx + mForces[1] * wFy + mForces[2] * wFz ) * jaci;
              double fy = -( mForces[1] * wFx + mForces[3] * wFy + mForces[4] * wFz ) * jaci;
              double fz = -( mForces[2] * wFx + mForces[4] * wFy + mForces[5] * wFz ) * jaci;

              if( 1 <= k && k <= Nz && ( fx != 0 || fy != 0 || fz != 0 ) )
              {
                GridPointSource* sourcePtr = new GridPointSource( mFreq, mT0, i, j, k, g,
                                                                  fx,
                                                                  fy, fz, mTimeDependence, mNcyc,
                                                                  mPar,
                                                                  mNpar, mIpar, mNipar );
                point_sources.push_back( sourcePtr );
              }
            }
          }
    }
  }
}

//-----------------------------------------------------------------------
void Source::exact_testmoments( int kx[3], int ky[3], int kz[3], double momex[3] )
{
  // Integrals over the domain of a polynomial of degree (kx,ky,kz) times the source
  if( !mIsMomentSource )
  {
    double x1, y1, z1;
    for( int c = 0 ; c < 3 ; c++ )
    {
      if( kx[c] == 0 )
        x1 = 1;
      else
        x1 = pow( mX0, kx[c] );
      if( ky[c] == 0 )
        y1 = 1;
      else
        y1 = pow( mY0, ky[c] );
      if( kz[c] == 0 )
        z1 = 1;
      else
        z1 = pow( mZ0, kz[c] );
      momex[c] = mForces[c] * x1 * y1 * z1;
    }
  }
  else
  {
    double x1, y1, z1, xp1, yp1, zp1;
    for( int c = 0 ; c < 3 ; c++ )
    {
      if( kx[c] == 0 )
        x1 = 1;
      else
        x1 = pow( mX0, kx[c] );
      if( kx[c] == 0 )
        xp1 = 0;
      else if( kx[c] == 1 )
        xp1 = -1;
      else
        xp1 = -kx[c] * pow( mX0, ( kx[c] - 1 ) );

      if( ky[c] == 0 )
        y1 = 1;
      else
        y1 = pow( mY0, ky[c] );
      if( ky[c] == 0 )
        yp1 = 0;
      else if( ky[c] == 1 )
        yp1 = -1;
      else
        yp1 = -ky[c] * pow( mY0, ( ky[c] - 1 ) );

      if( kz[c] == 0 )
        z1 = 1;
      else
        z1 = pow( mZ0, kz[c] );
      if( kz[c] == 0 )
        zp1 = 0;
      else if( kz[c] == 1 )
        zp1 = -1;
      else
        zp1 = -kz[c] * pow( mZ0, ( kz[c] - 1 ) );
      if( c == 0 )
        momex[c] = -( mForces[0] * xp1 * y1 * z1 + mForces[1] * x1 * yp1 * z1 + mForces[2] * x1 * y1 * zp1 );
      else if( c == 1 )
        momex[c] = -( mForces[1] * xp1 * y1 * z1 + mForces[3] * x1 * yp1 * z1 + mForces[4] * x1 * y1 * zp1 );
      else
        momex[c] = -( mForces[2] * xp1 * y1 * z1 + mForces[4] * x1 * yp1 * z1 + mForces[5] * x1 * y1 * zp1 );
    }
  }
}

//-----------------------------------------------------------------------
void Source::filter_timefunc( Filter* filter_ptr, double tstart, double dt, int nsteps )
{
  if( !m_is_filtered )
  {
    double (*timeFunc)( double f, double t, double* par, int npar, int* ipar, int nipar );
    switch( mTimeDependence )
    {
      case iRicker:
        timeFunc = RickerWavelet;
        break;
      case iGaussian:
        timeFunc = Gaussian;
        break;
      case iRamp:
        timeFunc = Ramp;
        break;
      case iTriangle:
        timeFunc = Triangle;
        break;
      case iSawtooth:
        timeFunc = Sawtooth;
        break;
      case iSmoothWave:
        timeFunc = SmoothWave;
        break;
      case iErf:
        timeFunc = Erf;
        break;
      case iVerySmoothBump:
        timeFunc = VerySmoothBump;
        break;
      case iC6SmoothBump:
        timeFunc = C6SmoothBump;
        break;
      case iRickerInt:
        timeFunc = RickerInt;
        break;
      case iBrune:
        timeFunc = Brune;
        break;
      case iBruneSmoothed:
        timeFunc = BruneSmoothed;
        break;
      case iDBrune:
        timeFunc = DBrune;
        break;
      case iGaussianWindow:
        timeFunc = GaussianWindow;
        break;
      case iLiu:
        timeFunc = Liu;
        break;
      case iDirac:
        timeFunc = Dirac;
        break;
      case iDiscrete:
        timeFunc = Discrete;
        break;
      default:
//        cout << "ERROR in Source::filter_timefunc, source type not recoginzed" << endl;
        throw GPException("ERROR in Source::filter_timefunc, source type not recoginzed");
    }

    // Convert to discrete representation
    double *discfunc = new double[nsteps];
    for( int k = 0 ; k < nsteps ; k++ )
      discfunc[k] = timeFunc( mFreq, tstart + k * dt - mT0, mPar, mNpar, mIpar, mNipar );
    mTimeDependence = iDiscrete;

// Filter the discretized function 
    filter_ptr->evaluate( nsteps, &discfunc[0], &discfunc[0] );

// Give the source time function a smooth start if this is a 2-pass (forward + backward) bandpass filter
    if( filter_ptr->get_passes() == 2 && filter_ptr->get_type() == bandPass )
    {
      double wghv, xi;
      int p0 = 3, p = 20; // First non-zero time level, and number of points in ramp;

      for( int i = 1 ; i <= p0 - 1 ; i++ )
      {
        discfunc[i - 1] = 0;
      }
      for( int i = p0 ; i <= p0 + p ; i++ )
      {
        wghv = 0;
        xi = ( i - p0 ) / ( (double) p );
        // polynomial P(xi), P(0) = 0, P(1)=1
        wghv = xi * xi * xi * xi * ( 35 - 84 * xi + 70 * xi * xi - 20 * xi * xi * xi );
        discfunc[i - 1] *= wghv;
      }
    }

    // Save discrete function
    mNipar = 1;
    mIpar = new int[mNipar];
    mIpar[0] = nsteps;

    mFreq = 1. / dt;
    delete[] mPar;
    mNpar = nsteps + 1;
    mPar = new double[mNpar];
    //      mPar[0] = tstart;
    mPar[0] = tstart - mT0;
    for( int i = 0 ; i < nsteps ; i++ )
      mPar[i + 1] = discfunc[i];
    delete[] discfunc;

    // Build the spline representation
    spline_interpolation();
    m_is_filtered = true;
  }
}

//-----------------------------------------------------------------------
int Source::spline_interpolation()
{
  // Assume mPar[1], to mPar[npts] contain the function
  // Assume mIpar[0] contains npts
  // Assume mFreq contains 1/dt, and mPar[0] is tstart.
  // Compute the six spline coefficients for each interval and return in mPar[1],to mPar[6*(npts-1)]
  if( mTimeDependence == iDiscrete )
  {
    int npts = mIpar[0];
    Qspline quinticspline( npts, &mPar[1], mPar[0], 1 / mFreq );
    double tstart = mPar[0];
    delete[] mPar;

    mNpar = 6 * ( npts - 1 ) + 1;
    mPar = new double[mNpar];
    mPar[0] = tstart;
    double* qsppt = quinticspline.get_polycof_ptr();
    for( int i = 0 ; i < 6 * ( npts - 1 ) ; i++ )
      mPar[i + 1] = qsppt[i];
    return 1;
  }
  else
    return 0;
}

//-----------------------------------------------------------------------
double Source::find_min_exponent() const
{
  // smallest number x, such that exp(x) does not cause underflow
  return -700.0;
}

//-----------------------------------------------------------------------
bool Source::is_inside( int i, int j, int k, int dims[6] ) const
                        {
  return dims[0] <= i && i <= dims[1] && dims[2] <= j && j <= dims[3]
      && dims[4] <= k && k <= dims[5];
}
