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

#ifndef SW4_MATERIALBLOCK_H
#define SW4_MATERIALBLOCK_H

#include "MaterialData.h"

class MaterialBlock : public MaterialData
{
public:
  MaterialBlock( double rho, double vs, double vp, double xmin, double xmax, double ymin,
                 double ymax, double zmin, double zmax, double qs=-1, double qp=-1,
                 double freq=1 );

  void set_material_properties( array<realT> &rho, array<realT> &cs,
                                array<realT> &cp, array<R1Tensor>& coord, realT zsurf );
  void set_gradients( double rhograd, double vsgrad, double vpgrad );
  void set_absoluteDepth( bool absDepth );
private:
  bool inside_block( double x, double y, double z );
  double m_rho, m_vp, m_vs, m_qp, m_qs, m_freq;
  double m_vpgrad, m_vsgrad, m_rhograd;
  double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
  double m_tol;
  bool m_absoluteDepth;
};

#endif
