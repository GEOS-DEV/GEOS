// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#include "MaterialBlock.h"

#include <iostream>

using namespace std;

//-----------------------------------------------------------------------
MaterialBlock::MaterialBlock( double rho, double vs, double vp, double xmin,
                              double xmax, double ymin, double ymax, double zmin, double zmax,
                              double qs, double qp, double freq )
{
  m_rho = rho;
  m_vp  = vp;
  m_vs  = vs;
  m_xmin = xmin;
  m_xmax = xmax;
  m_ymin = ymin;
  m_ymax = ymax;
  m_zmin = zmin;
  m_zmax = zmax;
  m_tol = 1e-5;
  m_vpgrad  = 0;
  m_vsgrad  = 0;
  m_rhograd = 0;
  m_qs = qs;
  m_qp = qp;
  m_freq = freq;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_absoluteDepth( bool absDepth )
{
  m_absoluteDepth = absDepth;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_gradients( double rhograd, double vsgrad, double vpgrad )
{
  m_rhograd = rhograd;
  m_vsgrad  = vsgrad;
  m_vpgrad  = vpgrad;
}

//-----------------------------------------------------------------------
bool MaterialBlock::inside_block( double x, double y, double z )
{
  return m_xmin-m_tol <= x && x <= m_xmax+m_tol && m_ymin-m_tol <= y &&
         y <= m_ymax+m_tol &&  m_zmin-m_tol <= z && z <= m_zmax+m_tol;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_material_properties( array<realT> & rho,
                                             array<realT>& cs, array<realT> & cp,
                                             array<R1Tensor>& coord, realT zsurf )
{
  for( localIndex ind = 0 ; ind < rho.size() ; ind++ )
    if(inside_block(coord[ind][0],coord[ind][1],coord[ind][2]))
    {
      if( m_rho != -1 )
        rho[ind] = m_rho + m_rhograd*(coord[ind][2]-zsurf);
      if( m_vs != -1 )
        cs[ind]  = m_vs + m_vsgrad*(coord[ind][2]-zsurf);
      if( m_vp != -1 )
        cp[ind]  = m_vp + m_vpgrad*(coord[ind][2]-zsurf);
    }
}
