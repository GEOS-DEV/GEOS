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

#ifndef SW4_MATERIALPFILE_H
#define SW4_MATERIALPFILE_H

#include <string>

#include "MaterialData.h"

using namespace std;

class MaterialPfile : public MaterialData
{
public:

  MaterialPfile( const std::string file,
                 const std::string directory,
                 const int nstenc,
                 const double vpminppm,
                 const double vsminppm,
                 const double rhominppm );

  ~MaterialPfile();

  virtual void set_material_properties( array<realT> &rho, array<realT> &cs,
                                        array<realT> &cp,  array<R1Tensor>& coord,
                                        realT zsurf );

  //  int get_material_pt( double x, double y, double z, double& rho, double&
  // cs, double& cp,
  //		       double& qs, double& qp );

  //  void getMinMaxBoundsZ(double& zmin, double& zmax);

protected:
  //  inline bool inside( double lat, double lon, double depth )
  //  {
  //    return m_latmin <= lat && lat <= m_latmax && m_lonmin <= lon && lon <=
  // m_lonmax
  //      && m_depthmin <= depth && depth <= m_depthmax;
  //  }

  inline bool inside_cart( double x, double y, double depth )
  {
    return m_xmin <= x && x <= m_xmax && m_ymin <= y && y <= m_ymax
           && m_depthmin <= depth && depth <= m_depthmax;
  }

  void read_pfile( );

  void sample_cart(double xs, double ys, double zs, double &vp,
                   double &vs, double &rho, double &qp, double &qs, bool debug );

  int m_nlat, m_nlon, m_nmaxdepth, m_nx, m_ny;
  int m_nstenc;
  double m_h, m_dlon, m_dlat;
  int     m_ksed, m_kmoho, m_k410, m_k660;
  double *m_x, *m_y;
// 3-dimensional arrays
  double* mZ, *mVp, *mVs, *mRho, *mQp, *mQs;

  double  m_vpmin, m_vsmin, m_rhomin;
  string m_model_file, m_model_dir, m_model_name;
  bool m_qf;

  double m_latmin, m_latmax, m_lonmin, m_lonmax, m_depthmin, m_depthmax;
  double m_xmin, m_xmax, m_ymin, m_ymax;
  //   bool m_coords_geographic;
};
#endif
