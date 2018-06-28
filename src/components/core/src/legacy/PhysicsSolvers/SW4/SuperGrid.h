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

#ifndef SW4SUPERGRID_H
#define SW4SUPERGRID_H
#include "Common/intrinsic_typedefs.h"

class SuperGrid
{
public:
  SuperGrid();
  void define_taper(bool left, realT leftStart, bool right, realT rightEnd,
                    realT width );
  realT dampingCoeff(realT x) const;
  realT stretching( realT x ) const;
  realT cornerTaper( realT x ) const;
  void   set_twilight( realT omega );
  void   print_parameters() const;

private:
  bool m_left, m_right;
  realT m_x0, m_x1, m_width, m_trans_width, m_const_width;
  realT m_epsL;
  realT sigma(realT xi) const;
  realT sigmaScale(realT x) const;
  realT linTaper(realT x) const;
};

#endif
