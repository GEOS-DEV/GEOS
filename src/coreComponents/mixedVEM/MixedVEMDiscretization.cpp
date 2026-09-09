/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MixedVEMDiscretization.cpp
 */

#include "mixedVEM/MixedVEMDiscretization.hpp"

namespace geos
{

using namespace dataRepository;

MixedVEMDiscretization::MixedVEMDiscretization( string const & name,
                                                Group * const parent )
  : Group( name, parent ),
  m_hybridization( 0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::hybridizationString(), &m_hybridization ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Solve the hybridized form of the mixed VEM system.\n"
                    "If 0, the indefinite saddle point system in the face tractions and the "
                    "element displacements is assembled and solved directly.\n"
                    "If 1, the stress space is broken elementwise, traction continuity is "
                    "restored by a Lagrange multiplier on the interior faces, and both element "
                    "unknowns are statically condensed. The global system is then the symmetric "
                    "positive definite interface problem H lambda = h, and the element stress and "
                    "displacement are recovered independently on each cell." );

  registerWrapper( viewKeyStruct::stabilizationLengthString(), &m_stabilizationLength ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Length h of the stabilization term, equation (15) of the reference.\n"
                    "If 0, the element diameter h_E is used, which is the choice of the paper.\n"
                    "If 1, the hydraulic radius |E| / |dE| is used. It is the only length whose "
                    "sum over the faces of h |f| is |E| for every shape and element type, so the "
                    "stabilization keeps its balance against the consistency term on flattened or "
                    "stretched cells, where h_E is the long diagonal of every face. It is more "
                    "accurate on such meshes and costs iterations, the stabilization being "
                    "smaller." );
}

MixedVEMDiscretization::CatalogInterface::CatalogType &
MixedVEMDiscretization::getCatalog()
{
  static MixedVEMDiscretization::CatalogInterface::CatalogType catalog;
  return catalog;
}

} // namespace geos
