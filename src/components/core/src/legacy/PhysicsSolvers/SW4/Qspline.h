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

#ifndef EW_QSPLINE_H
#define EW_QSPLINE_H

class Qspline
{
  int m_npts;
  double* m_polcof;
  double m_tmin, m_dt, m_dti;
public:
  Qspline( int npts, double* fun, double tmin, double dt, int bclow=1, int bchigh=1,
           double s1=0, double t1=0, double sn=0, double tn=0 );
  void Qsplineold( int npts, double* fun, double tmin, double dt );
  ~Qspline();
  double* get_polycof_ptr() { return m_polcof;}
  void evalf( double t, double& f );
  //   void evald( double t, double& f, double& f1, double& f2 );
  void evaldd( double t, double& f, double& f1, double& f2, double& f3, double& f4 );
};

#endif
