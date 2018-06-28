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

#ifndef GRID_POINT_SOURCE_H
#define GRID_POINT_SOURCE_H

#include <iostream>
#include "TimeDep.h"

class GridPointSource
{
  friend std ::ostream& operator<<(std::ostream& output, const GridPointSource& s);
public:

  GridPointSource(double frequency, double t0,
                  int i0, int j0, int k0, int g,
                  double Fx, double Fy, double Fz,
                  timeDep tDep, int ncyc, double* pars, int npar, int* ipars, int nipar );

  ~GridPointSource();

  int m_i0,m_j0,m_k0; // grid point index
  int m_grid;

  void getFxyz( double t, double* fxyz ) const;
  void getFxyztt( double t, double* fxyz ) const;
  void getFxyz_notime( double* fxyz ) const;

  double getTimeFunc(double t) const;
  double evalTimeFunc_t(double t) const;
  double evalTimeFunc_tt(double t) const;
  double evalTimeFunc_ttt(double t) const;
  double evalTimeFunc_tttt(double t) const;

  void limitFrequency(double max_freq);

  void print_info() const;

private:

  GridPointSource();

  void initializeTimeFunction();
  double mForces[3];
  double mFreq, mT0;

  timeDep mTimeDependence;
  double (*mTimeFunc)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_t)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_ttt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_om)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_omtt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tttt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tttom)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_ttomom)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tom)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_omom)(double f, double t,double* par, int npar, int* ipar, int nipar );

  double* mPar;
  int* mIpar;
  int  mNpar, mNipar;

  int mNcyc;
};

#endif
