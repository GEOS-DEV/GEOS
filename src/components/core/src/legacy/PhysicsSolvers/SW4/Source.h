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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SW4_SOURCE_H
#define SW4_SOURCE_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

class GridPointSource;
class Filter;

#include "TimeDep.h"

class Source
{
  friend std ::ostream& operator<<(std::ostream& output, const Source& s);
public:
  Source( double frequency, double t0, double x0, double y0, double z0,
          double Mxx, double Mxy, double Mxz, double Myy, double Myz, double Mzz,
          timeDep tDep, const char *name, bool topodepth, int ncyc=1,
          double* pars=NULL, int npars=0, int* ipars=NULL, int nipars=0 );

  Source( double frequency, double t0, double x0, double y0, double z0,
          double Fx, double Fy, double Fz,
          timeDep tDep, const char *name, bool topodepth, int ncyc=1,
          double* pars=NULL, int npars=0, int* ipars=NULL, int nipars=0 );
  ~Source();

  int m_i0, m_j0, m_k0;

  double getX0() const;
  double getY0() const;
  double getZ0() const;
  double getDepth() const;
  bool ignore() const {return mIgnore;}
  bool myPoint(){ return m_myPoint; }
  double getAmplitude() const;
  // Offset in time
  double getOffset() const;

  // Frequency
  double getFrequency() const;
  void setFrequency( double freq );
  timeDep getTfunc() const {return mTimeDependence;}
  void setMaxFrequency(double max_freq);

  // Type of source
  bool isMomentSource() const;

  double dt_to_resolve( int ppw ) const;
  int ppw_to_resolve( double dt ) const;

  const std::string& getName() const { return mName;};
  void limit_frequency( int ppw, double minvsoh );
  double compute_t0_increase( double t0_min ) const;
  void adjust_t0( double dt0 );

  void set_grid_point_sources4( double h, double xmin, double ymin, double zmin,
                                int Ni, int Nj, int Nz,
                                std::vector<GridPointSource*>& point_sources, int interior[6] ) const;

  void exact_testmoments( int kx[3], int ky[3], int kz[3], double momexact[3] );
  void getForces( double& fx, double& fy, double& fz ) const;
  void getMoments( double& mxx, double& mxy, double& mxz, double& myy, double& myz, double& mzz ) const;
  void setMoments( double mxx, double mxy, double mxz, double myy, double myz, double mzz );
  //      void printPointer(){std::cout << "Source pointer = "  << mPar <<
  // std::endl;}
  //  void perturb( double h, int comp );
  //  void set_derivative( int der );
  //  void set_noderivative( );
  //  void set_dirderivative( double dir[11] );
  //  Source* copy( std::string a_name );
  //  void set_parameters( double x[11] );

  //  void get_parameters( double x[11] ) const;
  void filter_timefunc( Filter* fi, double tstart, double dt, int nsteps );
  //  bool get_CorrectForMu(){return mShearModulusFactor;};
  //  void set_CorrectForMu(bool smf){mShearModulusFactor=smf;};

private:
  Source();
  void correct_Z_level();

  int spline_interpolation( );
  void getsourcewgh(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getsourcedwgh(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getsourcewghlow(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getsourcedwghlow(double ai, double wgh[6], double dwghda[6], double ddwghda[6] ) const;
  void getmetwgh( double alph, double wgh[8], double dwgh[8], double ddwgh[8], double dddwgh[8] ) const;
  void getmetdwgh( double alph, double wgh[8] ) const;
  void getmetwgh7( double ai, double wgh[7] ) const;
  void getmetdwgh7( double ai, double wgh[7] ) const;

  double find_min_exponent() const;
  bool is_inside( int i, int j, int k, int dims[6] ) const;

  std::string mName;
  std::vector<double> mForces;
  bool    mIsMomentSource;
  double  mFreq, mT0;
  bool    m_myPoint;
  bool    m_zRelativeToTopography;
  double  mX0,mY0,mZ0;
  double* mPar;
  int*    mIpar;
  int     mNpar, mNipar;
  int     mNcyc;
  timeDep mTimeDependence;
  bool    m_is_filtered;
  double  m_zTopo;
  bool    mIgnore;
};
//}

#endif
