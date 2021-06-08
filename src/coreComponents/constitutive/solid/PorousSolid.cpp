/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file PorousSolid.cpp
 */

#include "PorousSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "NullModel.hpp"
#include "porosity/BiotPorosity.hpp"
#include "porosity/PressurePorosity.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::PorousSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity >( name, parent )
{
}

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::~PorousSolid()
{}

// Register all CoupleSolid model types.
typedef PorousSolid< ElasticIsotropic, BiotPorosity > PoroElasticIsotropic;
typedef PorousSolid< ElasticTransverseIsotropic, BiotPorosity > PoroElasticTransverseIsotropic;
typedef PorousSolid< DruckerPrager, BiotPorosity > PoroDruckerPrager;
typedef PorousSolid< DruckerPragerExtended, BiotPorosity > PoroDruckerPragerExtended;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPragerExtended, string const &, Group * const )

}
} /* namespace geosx */
