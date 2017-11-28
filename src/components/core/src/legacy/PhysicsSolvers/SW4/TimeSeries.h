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
