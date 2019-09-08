/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file Tetrahedron.h
 */

#include "finiteElement/ElementLibrary/FiniteElementBase.h"

#ifndef TETRAHEDRON_H_
#define TETRAHEDRON_H_

namespace geosx
{
class TetrahedralElement : public FiniteElementBase
{
public:
  TetrahedralElement( BasisBase const & basis,
                      QuadratureBase const & quadrature,
                      const int num_zero_energy_modes );

  ~TetrahedralElement() override;

  static string CatalogName() { return "C3D4"; }

  void reinit( arrayView1d< R1Tensor const > const & X_global, arraySlice1d< localIndex const > const & mapped_support_points ) override;

};
}

#endif /* TETRAHEDRON_H_ */
