// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#include <mpi.h>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <time.h>

#include "TimeSeries.h"

using namespace std;

TimeSeries::TimeSeries( std::string fileName, std::string staName, receiverMode mode,
                        double x, double y, double z,
                        int i0, int j0, int k0, double xg, double yg, double zg,
                        int writeEvery, std::string path ):
  m_i0(i0),
  m_j0(j0),
  m_k0(k0),
  m_mode(mode),
  m_nComp(0),
  m_myPoint(true),
  m_fileName(fileName),
  m_staName(staName),
  mX(x),
  mY(y),
  mZ(z),
  mGPX(xg),
  mGPY(yg),
  mGPZ(zg),
  mWriteEvery(writeEvery),
  m_path(path),
  m_t0(),
  m_shift(0),
  m_dt(),
  mAllocatedSize(-1),
  mLastTimeStep(-1),
  mRecordedSol(),
  m_utc(),
  mQuietMode(false)
{
  if (!mQuietMode )
  {
    cout << "Receiver INFO for station " << m_fileName << ":" << endl <<
      "     initial location (x,y,z) = " << mX << " " << mY << " " << mZ << endl <<
      "     nearest grid point (x,y,z) = " << mGPX << " " << mGPY << " " << mGPZ <<
      " with indices (i,j,k)= " << m_i0 << " " << m_j0 << " " << m_k0 <<  endl;
  }

// get number of components from m_mode
  if (m_mode == Displacement || m_mode == Velocity)
    m_nComp=3;
  else if (m_mode == Div)
    m_nComp=1;
  else if (m_mode == Curl)
    m_nComp=3;
  else if (m_mode == Strains)
    m_nComp=6;

// allocate handles to solution array pointers
//   mRecordedSol = new double*[m_nComp];
//   for (int q=0; q<m_nComp; q++)
//      mRecordedSol[q] = static_cast<double*>(0);

// Set UTC as current date
  time_t tsec;
  time( &tsec );
  struct tm *utctime = gmtime( &tsec );
  m_utc[0] = utctime->tm_year+1900;
  m_utc[1] = utctime->tm_mon+1;
  m_utc[2] = utctime->tm_mday;
  m_utc[3] = utctime->tm_hour;
  m_utc[4] = utctime->tm_min;
  m_utc[5] = utctime->tm_sec;
  m_utc[6] = 0;  //milliseconds not given by 'time', not needed here.

} // end constructor

//--------------------------------------------------------------
TimeSeries::~TimeSeries()
{
// deallocate the recording arrays
//  if (mRecordedSol)
//  {
//    for (int q=0; q<m_nComp; q++)
//    {
//      if (mRecordedSol[q])
//	delete [] mRecordedSol[q];
//    }
//    delete [] mRecordedSol;
//  }
}

//--------------------------------------------------------------
void TimeSeries::allocateRecordingArrays( int numberOfTimeSteps,
                                          double startTime, double timeStep )
{
  if (numberOfTimeSteps > 0)
  {
    //  {
    //    mAllocatedSize = numberOfTimeSteps+1;
    mLastTimeStep = -1;
    //    for (int q=0; q<m_nComp; q++)
    //    {
    //      if (mRecordedSol[q]) delete [] mRecordedSol[q];
    //      mRecordedSol[q] = new double[mAllocatedSize];
    //    }
  }
  m_t0 = startTime;
  m_dt = timeStep;
}

//--------------------------------------------------------------
void TimeSeries::recordData( vector<realT> & u )
{
  if (!m_myPoint)
    return;

// better pass the right amount of data!
  if ( (int)u.size() != m_nComp )
  {
    cout << "SW4Solver TimeSeries::recordData: passing a vector of size= " << u.size() << " but nComp= " <<  m_nComp << endl;
    return;
  }

// ---------------------------------------------------------------
// This routine only knows how to push the nComp doubles on the array stack.
// The calling routine need to figure out what needs to be saved
// and do any necessary pre-calculations
// ---------------------------------------------------------------
  mLastTimeStep++;
  for( int q=0 ; q <m_nComp ; q++ )
    mRecordedSol.push_back(u[q]);
  //   if (mLastTimeStep < mAllocatedSize)
  //   {
  //      for (int q=0; q<m_nComp; q++)
  //	 mRecordedSol[q][mLastTimeStep] = u[q];
  //   }
  //   else
  //   {
  //      cout << "SW4Solver, TimeSeries: Ran out of recording space for the
  // receiver station at (i,j,k) =  "
  //	   << m_i0 << " " << m_j0 << " " << m_k0 << endl;
  //      return;
  //   }
  if (mWriteEvery > 0 && mLastTimeStep > 0 && mLastTimeStep % mWriteEvery == 0)
    writeFile();
}


//----------------------------------------------------------------------
void TimeSeries::writeFile( string suffix )
{
  if (!m_myPoint)
    return;

// ------------------------------------------------------------------
// We should add an argument to this function that describes how the
// header and filename should be constructed
// ------------------------------------------------------------------

  stringstream filePrefix;

//building the file name...
  if( m_path != "." && m_path != "./" )
  {
    filePrefix << m_path;
    if( m_path[m_path.size()-1] != '/' )
      filePrefix << "/";
  }
  if( suffix == "" )
    filePrefix << m_fileName << ".";
  else
    filePrefix << m_fileName << suffix.c_str() << ".";

// Write out displacement components (ux, uy, uz)
  filePrefix << "txt";
  if ( !mQuietMode )
    cout << "Writing ASCII USGS file, "
         << "of size " << mLastTimeStep+1 << ": "
         << filePrefix.str() << endl;
  write_usgs_format( filePrefix.str() );
}

//-----------------------------------------------------------------------
void TimeSeries::write_usgs_format(string a_fileName)
{
  //  string mname[] =
  // {"Zero","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
  FILE *fd=fopen(a_fileName.c_str(),"w");
  //  double lat, lon;
  //  double x, y, z;

  if( fd == NULL )
    cout << "ERROR: opening USGS file " << a_fileName << " for writing" <<  endl;
  else
  {

    //      cout << "IN write_usgs " << " mt0 = " << m_t0 << " mshift = " <<
    // m_shift << endl;
// frequency resolution
//    double freq_limit=-999;
//    if (a_ew->m_prefilter_sources)
//      freq_limit = a_ew->m_fc;
//    else if (a_ew->m_limit_frequency)
//      freq_limit = a_ew->m_frequency_limit;

// write the header
    fprintf(fd, "# Author: SW4\n");
    fprintf(fd, "# Scenario: %s\n", "test" /*a_ew->m_scenario.c_str()*/);
//   if( m_utc_set )
// AP: micro-second field is padded from left with 0, i.e., 1 micro sec gets
// written as 001, second is also padded by a zero, if needed
    fprintf(fd, "# Date: UTC  %02i/%02i/%i:%i:%i:%02i.%.3i\n", m_utc[1], m_utc[2], m_utc[0], m_utc[3],
            m_utc[4], m_utc[5], m_utc[6] );
    //   else
    //      fprintf(fd, "# Date: %i-%s-%i\n", mEventDay,
    // mname[mEventMonth].c_str(), mEventYear);

    fprintf(fd, "# Bandwith (Hz): %e\n", 1.234 /*freq_limit*/);
    fprintf(fd, "# Station: %s\n", m_fileName.c_str() /*mStationName.c_str()*/ );
    //      fprintf(fd, "# Target location (WGS84 longitude, latitude) (deg): %e
    // %e\n", m_rec_lon, m_rec_lat);
    //      fprintf(fd, "# Actual location (WGS84 longitude, latitude) (deg): %e
    // %e\n", m_rec_gp_lon, m_rec_gp_lat);
    fprintf(fd, "# Target location (x,y,z) : %e %e %e \n", mX, mY, mZ);
    fprintf(fd, "# Actual location (x,y,z) : %e %e %e \n", mGPX, mGPY, mGPZ );
// distance in horizontal plane
    fprintf(fd, "# Distance from target to actual location (m): %e\n", sqrt( (mX-mGPX)*(mX-mGPX)+(mY-mGPY)*(mY-mGPY) ) );
    fprintf(fd, "# nColumns: %i\n", m_nComp+1);

    fprintf(fd, "# Column 1: Time (s)\n");
    if (m_mode == Displacement )
    {
      fprintf(fd, "# Column 2: X displacement (m)\n");
      fprintf(fd, "# Column 3: Y displacement (m)\n");
      fprintf(fd, "# Column 4: Z displacement (m)\n");
    }
    else if( m_mode == Velocity )
    {
      fprintf(fd, "# Column 2: X velocity (m/s)\n");
      fprintf(fd, "# Column 3: Y velocity (m/s)\n");
      fprintf(fd, "# Column 4: Z velocity (m/s)\n");
    }
    else if( m_mode == Div )
    {
      fprintf(fd, "# Column 2: divergence of displacement ()\n");
    }
    else if( m_mode == Curl )
    {
      fprintf(fd, "# Column 2: curl of displacement, component 1 ()\n");
      fprintf(fd, "# Column 3: curl of displacement, component 2 ()\n");
      fprintf(fd, "# Column 4: curl of displacement, component 3 ()\n");
    }
    else if( m_mode == Strains )
    {
      fprintf(fd, "# Column 2: xx strain component ()\n");
      fprintf(fd, "# Column 3: yy strain component ()\n");
      fprintf(fd, "# Column 4: zz strain component ()\n");
      fprintf(fd, "# Column 5: xy strain component ()\n");
      fprintf(fd, "# Column 6: xz strain component ()\n");
      fprintf(fd, "# Column 7: yz strain component ()\n");
    }

// write the data
    for( int i = 0 ; i <= mLastTimeStep ; i++ )
    {
      fprintf(fd, "%e", m_shift + i*m_dt);
      for (int q=0 ; q<m_nComp ; q++)
        fprintf(fd, " %20.12g", mRecordedSol[q+m_nComp*i]);
      fprintf(fd, "\n");
    }
    fclose(fd);
  }
}
