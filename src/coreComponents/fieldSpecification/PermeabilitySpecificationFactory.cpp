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
 * @file PermeabilitySpecificationFactory.cpp
 */

#include "PermeabilitySpecificationFactory.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "PermeabilitySpecification.hpp"
#include "FieldSpecificationABC.hpp"
#include "FieldSpecificationBase.hpp"

namespace geos
{

/**
 * @class PermeabilitySpecificationFactory
 */

void PermeabilitySpecificationFactory::generate( FieldSpecificationABC const & spec,
                                                 dataRepository::Group & manager ) const
{
  auto ps = dynamic_cast< PermeabilitySpecification const * >( &spec ); 

  stdArray< string, 3 > suffixes = {{ "_x", "_y", "_z" }};

  R1Tensor scales = ps->getScales();

  for ( string const & regionName : ps->getRegionNames() )
  {
    string const objectPath = "ElementRegions/" + regionName;

    for ( integer comp = 0; comp < 3; ++comp )
    {
      string const childName = ps->getName() + "_" + regionName + suffixes[ comp ];

      FieldSpecificationBase & fs = manager.registerGroup< FieldSpecificationBase >( childName );
      fs.setFieldName( ps->getFieldName() );
      fs.setObjectPath( objectPath );
      fs.setScale( scales[ comp ] );
      fs.initialCondition( true );
      fs.setComponent( comp );

      for ( auto const & setName : ps->getSetNames() )
      {
        fs.addSetName( setName );
      }

      if ( !ps->getFunctionName().empty() )
      {
        fs.setFunctionName( ps->getFunctionName() );
      }

    }
  }

}

}