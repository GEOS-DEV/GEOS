/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
