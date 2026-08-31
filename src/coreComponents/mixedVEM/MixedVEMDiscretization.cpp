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
}

MixedVEMDiscretization::CatalogInterface::CatalogType &
MixedVEMDiscretization::getCatalog()
{
  static MixedVEMDiscretization::CatalogInterface::CatalogType catalog;
  return catalog;
}

} // namespace geos
