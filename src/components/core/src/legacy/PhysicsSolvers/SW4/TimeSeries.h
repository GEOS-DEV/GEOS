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

#ifndef SW4_TIMESERIES_H
#define SW4_TIMESERIES_H

#include <vector>
#include <string>
#include "Common/typedefs.h"

using namespace std;

class TimeSeries
{

public:
  enum receiverMode {Displacement, Div, Curl, Strains, Velocity };

  TimeSeries( std::string fileName, std::string staName, receiverMode mode,
              double x, double y, double z, int i0, int j0, int k0,
              double xg, double yg, double zg, int writeEvery, std::string path="." );
  ~TimeSeries();

  void allocateRecordingArrays( int numberOfTimeSteps, double startTime, double timeStep );

  void recordData( vector<realT> & u);

  void writeFile( string suffix="" );

  int getNsteps() const {return mLastTimeStep+1;}

  bool myPoint(){ return m_myPoint; }

  receiverMode getMode(){ return m_mode; }

  double getX() const {return mX;}
  double getY() const {return mY;}
  double getZ() const {return mZ;}

  std::string getStationName(){return m_staName;}

// for simplicity, make the grid point location public
  int m_i0;
  int m_j0;
  int m_k0;

private:
  TimeSeries();

  void write_usgs_format( string a_fileName);

  //   void getwgh( double ai, double wgh[6], double dwgh[6], double ddwgh[6] );
  //   void getwgh5( double ai, double wgh[6], double dwgh[6], double ddwgh[6]
  // );

  receiverMode m_mode;
  int m_nComp;

  bool m_myPoint;  // set to true if this processor writes to the arrays

  std::string m_fileName, m_staName;

  double mX, mY, mZ, mGPX, mGPY, mGPZ;  // original and actual location

  int mWriteEvery;

  string m_path;

// start time, shift, and time step
  double m_t0, m_shift, m_dt;

// size of recording arrays
  int mAllocatedSize;

// index of last recorded element
  int mLastTimeStep;

// recording arrays
  std::vector<double> mRecordedSol;

//   double** mRecordedSol;

// UTC time for start of seismogram,
//     m_t0  =  m_utc - 'utc reference time',
//           where 'utc reference time' corresponds to simulation time zero.
//     the time series values are thus given by simulation times t_k = m_t0 +
// m_shift + k*m_dt, k=0,1,..,mLastTimeStep
  int m_utc[7];

  bool mQuietMode;
};

#endif
