/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
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
