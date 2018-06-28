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

#include <cstring>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "MaterialPfile.h"

using namespace std;

//-----------------------------------------------------------------------
MaterialPfile::MaterialPfile( const std::string file,
                              const std::string directory,
                              const int nstenc,
                              const double vpmin,
                              const double vsmin,
                              const double rhomin )
  :
  m_nlat(0),
  m_nlon(0),
  m_nmaxdepth(0),
  m_nx(0),
  m_ny(0),
  m_nstenc(nstenc),
  m_h(0.0),
  m_dlon(0.0),
  m_dlat(0.0),
  m_ksed(0),
  m_kmoho(0),
  m_k410(0),
  m_k660(0),
  m_x(nullptr),
  m_y(nullptr),
  mZ(nullptr),
  mVp(nullptr),
  mVs(nullptr),
  mRho(nullptr),
  mQp(nullptr),
  mQs(nullptr),
  m_vpmin(vpmin),
  m_vsmin(vsmin),
  m_rhomin(rhomin),
  m_model_file(file),
  m_model_dir(directory),
  m_model_name()
{
  read_pfile();
}

//-----------------------------------------------------------------------
void MaterialPfile::set_material_properties( array<realT>& rho,
                                             array<realT>& cs,
                                             array<realT>& cp,
                                             array<R1Tensor>& coord,
                                             realT zsurf )
{
  int outside = 0; int material = 0;
  double vp,vs,density,qup,qus;
  for( localIndex ind = 0 ; ind < rho.size() ; ind++ )
  {
    if( inside_cart(coord[ind][0],coord[ind][1],coord[ind][2]) )
    {
      sample_cart( coord[ind][0], coord[ind][1], coord[ind][2]-zsurf,
                   vp, vs, density, qup, qus, false );
      rho[ind] = density;
      cs[ind] = vs;
      cp[ind] = vp;
      material++;
    }
    else
      outside++;
  }
  //   int outsideSum, materialSum;
  //   MPI_Reduce(&outside, &outsideSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD
  // );
  //   MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0,
  // MPI_COMM_WORLD );
  //   if (mEW->proc_zero())
  //      cout << "outside = " << outsideSum << ", " << "material = " <<
  // materialSum << endl;
}

//-----------------------------------------------------------------------
void MaterialPfile::read_pfile( )
{
  //   int kk, m;

  //   m_qf = false;

  //   m_model_dir  = ppdir;
  //   m_model_file = ppfile;
  string ppmfile = m_model_dir + "/" + m_model_file;

  //   m_vpmin   = vpmin;
  //   m_vsmin   = vsmin;
  //   m_rhomin  = rhomin;

  //   int myRank;
  //   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

// Open file
  FILE* fd=fopen(ppmfile.c_str(), "r" );
  if( fd == NULL )
  {
    cout << "Unable to open the pfile input file: '" << ppmfile << "'" << endl;
    //     if (myRank == 0) cerr << "Unable to open the pfile input file: '" <<
    // ppmfile << "'" << endl;
    //     MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int bufsize = 1024;
  int nread;
  char* buf = new char[bufsize];

  // Read pfile header

  // Line 1
  if( fgets(buf,bufsize,fd) == NULL )
    cout << "Error line 1 in pfile header not found\n";
  string tok0 = strtok( buf, " \t" );
// strip off any white space
  size_t nWhite = tok0.find_first_of(" \t\n");
  m_model_name = tok0.substr(0,nWhite);

  // Line 2
  if( fgets(buf,bufsize,fd) == NULL )
    cout << "Error line 2 in pfile header not found\n";
  nread = sscanf(buf,"%le", &m_h);
  if( nread!=1 )
    cout <<  "Error reading 2nd line of header, nread= " << nread << " but expected 1\n";

  // Line 3
  if( fgets(buf,bufsize,fd) == NULL )
    cout << "Error line 3 in pfile header not found\n";
  nread = sscanf(buf,"%i %le %le", &m_nlat, &m_latmin, &m_latmax);
  if( nread!=3 )
    cout << "Error reading 3rd line of header, nread= " << nread << " but expected 3\n";

  // Line 4
  if( fgets(buf,bufsize,fd) == NULL )
    cout << "Error line 4 in pfile header not found\n";
  nread = sscanf(buf,"%i %le %le", &m_nlon, &m_lonmin, &m_lonmax);
  if( nread!=3 )
    cout << "Error reading 4th line of header, nread= " << nread << " but expected 3\n";

  // Line 5
  if( fgets(buf,bufsize,fd) == NULL )
    cout << "Error line 5 in pfile header not found\n";
  nread = sscanf(buf,"%i %le %le", &m_nmaxdepth, &m_depthmin, &m_depthmax);
  if( nread!=3 )
    cout << "Error reading 5th line of header, nread= " << nread << " but expected 3\n";

  // Line 6
  if( fgets(buf,bufsize,fd) == NULL )
    cout << "Error line 6 in pfile header not found\n";
  nread = sscanf(buf,"%i %i %i %i", &m_ksed, &m_kmoho, &m_k410, &m_k660);
  if( nread!=4 )
    cout << "Error reading 6th line of header, nread= " << nread << " but expected 4\n";

  // Line 7
  if( fgets(buf,bufsize,fd) == NULL )
    cout << "Error line 7 in pfile header not found\n";
  char* tok = strtok(buf," \t");
  if( tok == NULL )
    cout << "Error on line 7 in pfile header, no Q-available flag\n";
  string cqf0 = tok;
// strip off any white space
  nWhite = cqf0.find_first_of(" \t\n");
  string cqf = cqf0.substr(0,nWhite);

// test
//   printf("Q-flag string '%s'\n", cqf.c_str());
  m_qf = ( cqf == "T") || (cqf == "t") || (cqf == ".TRUE.") || (cqf == ".true.");

// done reading the header

// make sure the stencil width does not exceed the number of grid points
  if (m_nstenc > m_nlat || m_nstenc > m_nlon)
  {
    m_nstenc = min(m_nlat, m_nlon);
    cout << "Warning: pfile: stencil width reduced to " << m_nstenc << endl;
  }

  m_nx = m_nlat;
  m_ny = m_nlon;
  m_xmin = m_latmin;
  m_xmax = m_latmax;
  m_ymin = m_lonmin;
  m_ymax = m_lonmax;

  //   cout << "Pfile model name (string): '" << m_model_name << "'" << endl;
  //   cout << "Step size in x and y: " << m_h << endl;
  //   cout << "Number of x-direction points: " << m_nx << endl;
  //   cout << "Min x: " << m_xmin << " Max x: " << m_xmax << endl;
  //   cout << "Number of y-direction points: " << m_ny << endl;
  //   cout << "Min y: " << m_ymin << " Max y: " << m_ymax << endl;
  //   cout << "Number of depth points: " << m_nmaxdepth << endl;
  //   cout << "Min depth: " << m_depthmin << " Max depth: " << m_depthmax <<
  // endl;
  //   cout << "Optional indices: Sediment: " << m_ksed << " MoHo: " << m_kmoho
  // << " 410: " << m_k410 << " 660: " << m_k660 << endl;
  //   cout << "Attenuation Q-factors available: " << (m_qf? "yes":"no") <<
  // endl;


  // Allocate arrays
  m_x = new double[m_nx];
  m_y = new double[m_ny];

// new 3-D Sarrays
//   mZ.define(m_nlon, m_nlat, m_nmaxdepth);
//   mVp.define(m_nlon, m_nlat, m_nmaxdepth);
//   mVs.define(m_nlon, m_nlat, m_nmaxdepth);
//   mRho.define(m_nlon, m_nlat, m_nmaxdepth);

  size_t npts = (static_cast<size_t>(m_nlon))*m_nlat*m_nmaxdepth;
  mZ = new double[npts];
  mVp = new double[npts];
  mVs = new double[npts];
  mRho = new double[npts];

  if(m_qf)
  {
// new 3-D Sarrays
    mQp = new double[npts];
    mQs = new double[npts];
  }
  //   else
  //   {
  //      cout << "ppmod: NOT allocating arrays for Qp and Qs\n";
  //   }

  int kk, ndepth, line=7;
  double zc, vp, vs, rho, qp, qs;
// note: nx = nlat, ny = nlon
  for(int jy=0 ; jy < m_ny ; jy++ )
    for(int ix=0 ; ix < m_nx ; ix++ )
    {
      if( fgets(buf,bufsize,fd) == NULL )
        cout << "Error in pfile profile header at coordinate " << ix << " " << jy << "\n";
      nread = sscanf(buf,"%le %le %i", &m_x[ix], &m_y[jy], &ndepth);
      if( !(nread==3) )
        cout << "Error reading 1st line of profile at " << ix << " " << jy
             << " nread= " << nread << " but expected 3\n";
      line++;
// fundamental sanity checks
      if (!(m_y[jy] <= m_ymax && m_y[jy] >= m_ymin &&
            m_x[ix] <= m_xmax && m_x[ix] >= m_xmin) )
      {
        cout << "Error reading pfile: x profile #" << ix+1 << " y profile #" << jy+1 <<
          ": x=" << m_x[ix] << " or " << m_y[jy] << " out of bounds! min(x)=" << m_xmin << " max(x)="
             << m_xmax << ", min(y)=" << m_ymin << ", max(y)=" << m_ymax << "\n";
      }
// sanity check 2
      if (ndepth != m_nmaxdepth )
      {
        //	       if (myRank == 0)
        //	       {
        cerr << "pfile reader error, ppmod file line=" << line << endl;
        cerr << "read ndepth=" << ndepth << " which is different from header nmaxdepth="
             << m_nmaxdepth << endl;
        //	       }
      }
// check the range
      double y = m_ymin + jy*m_h;
      double x = m_xmin + ix*m_h;
      if (fabs(y - m_y[jy]) + fabs(x - m_x[ix]) > 0.1*m_h)
      {
        //	       if (myRank == 0)
        //	       {
        cerr << "pfile reader error, ppmod file line=" << line << endl;
        cerr << "read x[" << ix << "]=" << m_x[ix] << " but expected x=" << x << endl;
        cerr << "read y[" << jy << "]=" << m_y[jy] << " but expected y=" << y << endl;
        cerr << "CHECK THE PPMOD FILE." << endl;
        cerr << "DEPTH PROFILES SHOULD BE ORDERED SUCH THAT X VARIES THE FASTEST" << endl;
        //	       }
        //	       MPI_Abort(MPI_COMM_WORLD, 1);
      }

// Read depth profile
      for(int k=0 ; k < m_nmaxdepth ; k++ )
      {
        if( fgets( buf, bufsize, fd ) == NULL )
          cout << "Error in pfile profile at coordinate " << ix << " " <<
            jy << " " << k << "\n";
        if (m_qf)
        {
          nread=sscanf(buf, "%i %le %le %le %le %le %le", &kk, &zc, &vp, &vs, &rho, &qp, &qs);
          if( nread !=7 )
            cout << "Error reading pfile at " << ix << " " << jy << " " << k
                 << " nread= " << nread << " but expected 7\n";
        }
        else
        {
          nread=sscanf(buf, "%i %le %le %le %le", &kk, &zc, &vp, &vs, &rho);
          if( nread != 5 )
            cout << "Error reading pfile at " << ix << " " << jy << " " << k
                 << " nread= " << nread << " but expected 5\n";
        }
        size_t ind = ix + m_nlon*jy + m_nlon*m_nlat*k;
        mZ[ind] = zc;
        mVp[ind] = vp;
        mVs[ind] = vs;
        mRho[ind] = rho;
        //         mZ(ix+1,jy+1,k+1) = zc;
        //         mVp(ix+1,jy+1,k+1) = vp;
        //         mVs(ix+1,jy+1,k+1) = vs;
        //         mRho(ix+1,jy+1,k+1) = rho;

        if (m_qf)
        {
          mQp[ind] = qp;
          mQs[ind] = qs;
          //		 mQp(ix+1,jy+1,k+1) = qp;
          //		 mQs(ix+1,jy+1,k+1) = qs;
        }
        line++;
      }
    }
  fclose(fd);

  //   if (myRank == 0)
  //   {
  //      cout << "******* Done reading Pfile **********" << endl << endl;
  //   }
  delete[] buf;
}

//-----------------------------------------------------------------------
void MaterialPfile::sample_cart( double xs, double ys, double zs, double &vp,
                                 double &vs, double &rho, double &qp, double &qs, bool debug )
//--------------------------------------------------------------------------
// return material properties (vp, vs, rho) at point (xs, ys, zs)
//--------------------------------------------------------------------------
{
// tmp
//  if (debug) printf("DEBUG::sample_cart: xs=%e, ys=%e, zs=%e, m_h=%e\n", xs,
// ys, zs, m_h);

  //  Check if xs and ys are out of range
  if ( xs < m_xmin )
  {
    cerr << "MaterialPfile::sample xs out of range (min): " <<  xs << ", " <<  m_xmin << endl;
    //     MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if ( xs > m_xmax )
  {
    cerr << "MaterialPfile::sample xs out of range (max): " << xs << ", " << m_xmax << endl;
    //     MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if ( ys < m_ymin)
  {
    cerr << "MaterialPfile::sample ys out of range (min): " << ys << ", " <<  m_ymin << endl;
    //     MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if( ys > m_ymax)
  {
    cerr << "MaterialPfile::sample ys out of range (max): " << ys << ", " << m_ymax << endl;
    //     MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int ii, jj;
  if( m_nstenc % 2 == 1 )
  {
    // odd number of points in stencil
    int s = (m_nstenc-1)/2;
    ii  = static_cast<int>( floor( (xs-m_xmin)/m_h +0.5) )-s;
    jj  = static_cast<int>( floor( (ys-m_ymin)/m_h +0.5) )-s;
  }
  else
  {
    // even number of points in stencil
    int s = m_nstenc/2;
    ii  = static_cast<int>( floor( (xs-m_xmin)/m_h ) )-s+1;
    jj  = static_cast<int>( floor( (ys-m_ymin)/m_h ) )-s+1;
  }

// make sure we stay within array boundaries
  if( ii < 0 )
    ii = 0;
  if( jj < 0 )
    jj = 0;

  int ii2 = ii + m_nstenc-1;
  int jj2 = jj + m_nstenc-1;

// m_x = new double[m_nx], m_nx=m_nlat
  if( ii2 >= m_nx )
  {
    ii2 = m_nx-1;
    ii = ii2 -(m_nstenc-1);
  }

// m_y = new double[m_ny], m_ny=m_nlon
  if( jj2 >= m_ny )
  {
    jj2 = m_ny-1;
    jj = jj2 - (m_nstenc-1);
  }

  double w=0;
  vp=vs=rho=qp=qs=0;
  double appm  = 0.5*m_nstenc*m_h/sqrt(-log(1e-6));
  double appmi2 = 1.0/(appm*appm);

  for( int j1 = jj ; j1 <= jj2 ; j1++ )
    for( int i1 = ii ; i1 <= ii2 ; i1++ )
    {
      double wgh = exp(-( (xs-m_x[i1])*(xs-m_x[i1])
                          +(ys-m_y[j1])*(ys-m_y[j1]) )*appmi2 );
      w += wgh;

// depth index
      int kk;
      for( kk=1 ; kk < m_nmaxdepth ; kk++ )// AP changed from kk <= m_nmaxdepth
      {
        size_t ind = i1 + m_nx*j1 + m_nx*m_ny*(kk-1);
        //	     if (mZ(i1+1,j1+1,kk) > zs) break;
        if (mZ[ind] > zs)
          break;
      }
// at this point we should have mZ(kk-1) <= zs < mZ(kk), kk <= m_nmaxdepth

      int k1 = kk-1;
// now we should have mZ(k1) <= zs < mZ(k1+1)
      size_t ind  = i1 + m_nx*j1 + m_nx*m_ny*(k1-1);
      size_t indp = i1 + m_nx*j1 + m_nx*m_ny*(k1);
// linear interpolation factor
      double factor = (zs-mZ[ind])/(mZ[indp]-mZ[ind]);

// new style
      vp  += (mVp[ind]  + factor*(mVp[indp]-  mVp[ind]) )*wgh;
      vs  += (mVs[ind]  + factor*(mVs[indp]-  mVs[ind]) )*wgh;
      rho += (mRho[ind] + factor*(mRho[indp]-mRho[ind]) )*wgh;
// tmp
      // if (debug) printf("DEBUG: i1+1=%i, j1+1=%i, k1=%i, vp=%e, wgh=%e\n",
      // i1+1, j1+1, k1,
      //        mVp[ind] + factor*(mVp[indp]-mVp[ind]), wgh);

      if( m_qf )
      {
        // qp += (m_qp[m+k1] + factor*(m_qp[m+k1+1]-m_qp[m+k1]))*wgh;
        // qs += (m_qs[m+k1] + factor*(m_qs[m+k1+1]-m_qs[m+k1]))*wgh;
        qp += (mQp[ind] + factor*(mQp[indp]-mQp[ind]) )*wgh;
        qs += (mQs[ind] + factor*(mQs[indp]-mQs[ind]) )*wgh;
      }
    }

// Normalize
  double iw = 0.0;
  if (w != 0.)
    iw = 1.0/w;
  else
  {
    printf("Error MaterialPfile::sample_cart: weight w = 0 at x=%e, y=%e, depth=%e\n", xs, ys, zs);
// tmp
    printf("ii=%i, ii2=%i, jj=%i, jj2=%i, x[ii]=%e, y[jj]=%e\n", ii, ii2, jj, jj2, m_x[ii], m_y[jj]);
    double dist2    = (xs-m_x[ii])*(xs-m_x[ii]) +(ys-m_y[jj])*(ys-m_y[jj]);
    double exponent = -dist2*appmi2;
    double wgh      = exp(exponent);
    printf("dist2=%e, appm=%e, appmi2=%e, exponent=%e, wgh=%e\n", dist2, appm, appmi2, exponent, wgh);
    //     MPI_Abort(MPI_COMM_WORLD, 1);
  }
  vp  *= iw;
  vs  *= iw;
  rho *= iw;
  if( m_qf )
  {
    qp *= iw;
    qs *= iw;
  }
}

//-----------------------------------------------------------------------
MaterialPfile::~MaterialPfile()
{
  delete[] m_x;
  delete[] m_y;
  delete[] mRho;
  delete[] mVp;
  delete[] mVs;
  delete[] mZ;
  if(m_qf)
  {
    delete[] mQs;
    delete[] mQp;
  }
}
