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

#include "CopyField.hpp"
#include "managers/DomainPartition.hpp"
#include "FieldSpecificationBase.hpp"

namespace geosx
{

using namespace dataRepository;

CopyField::CopyField( std::string const & name,
                        Group * const parent ):
  TaskBase( name, parent )
{
  registerWrapper(viewKeyStruct::fromString, &m_from, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Name of the field to copy.");

  registerWrapper(viewKeyStruct::toString, &m_to, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Name of the field that will be created and filled.");

  registerWrapper(viewKeyStruct::targetRegionsString, &m_targetRegions, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Regions on which the copy will be done.");
}

CopyField::~CopyField()
{
}


void CopyField::Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                         real64 const GEOSX_UNUSED_PARAM( dt ),
                         integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                         integer const GEOSX_UNUSED_PARAM( eventCounter ),
                         real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                         Group * domain )
{
  GEOSX_LOG_RANK( "Copying " + m_from + " to " + m_to );
  DomainPartition * domainCast = domain->group_cast<DomainPartition*>(domain);
  MeshBody * meshBody = domainCast->getMeshBody(0);
  MeshLevel * meshLevel = meshBody->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  for( auto regionName : m_targetRegions )
  {
    ElementRegionBase * curRegion = elemManager->GetRegion(regionName);
    GEOSX_ERROR_IF( curRegion == nullptr, "Region " + regionName + " not found" );
    curRegion->forElementSubRegions([&]( auto * const elementSubRegion) -> void
    {
      dataRepository::WrapperBase * vwFrom = elementSubRegion->getWrapperBase( m_from );
      GEOSX_ERROR_IF( vwFrom == nullptr, "Field " + m_from + " not found on " + elementSubRegion->getName() );
      std::type_index typeIndex = std::type_index( vwFrom->get_typeid());
      rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ),
                                      [&]( auto type ) -> void
      {
        using fieldType = decltype(type);
        dataRepository::Wrapper< fieldType > *  to = elementSubRegion->template registerWrapper< fieldType >(m_to);
        dataRepository::Wrapper< fieldType > & from = dataRepository::Wrapper< fieldType >::cast( *vwFrom );
        fieldType & fromField = from.reference();
        fieldType & toField = to->reference();
        toField = fromField;
      });
    });
  }

}

REGISTER_CATALOG_ENTRY( TaskBase, CopyField, std::string const &, Group * const )

} /* namespace */
