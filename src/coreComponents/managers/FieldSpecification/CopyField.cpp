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
                        ManagedGroup * const parent ):
  TaskBase( name, parent )
{

  RegisterViewWrapper(viewKeyStruct::fromString, &m_from, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Name of the field to copy.");

  RegisterViewWrapper(viewKeyStruct::toString, &m_to, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Name of the field that will be created and filled.");

  RegisterViewWrapper(viewKeyStruct::targetRegionsString, &m_targetRegions, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "Regions on which the copy will be done.");
}

CopyField::~CopyField()
{
}


void CopyField::Execute( real64 const time_n,
                         real64 const dt,
                         integer const cycleNumber,
                         integer const eventCounter,
                         real64 const eventProgress,
                         ManagedGroup * domain )
{
  GEOS_LOG_RANK( "Copying " + m_from + " to " + m_to );
  DomainPartition * domainCast = domain->group_cast<DomainPartition*>(domain);
  MeshBody * meshBody = domainCast->getMeshBody(0);
  MeshLevel * meshLevel = meshBody->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  for( auto regionName : m_targetRegions )
  {
    ElementRegion * curRegion = elemManager->GetRegion(regionName);
    GEOS_ERROR_IF( curRegion == nullptr, "Region " + regionName + " not found" );
    curRegion->forElementSubRegions([&]( auto * const elementSubRegion) -> void
    {
      dataRepository::ViewWrapperBase * vwFrom = elementSubRegion->getWrapperBase( m_from );
      GEOS_ERROR_IF( vwFrom == nullptr, "Field " + m_from + " not found on " + elementSubRegion->getName() );
      std::type_index typeIndex = std::type_index( vwFrom->get_typeid());
      rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ),
                                      [&]( auto type ) -> void
      {
        using fieldType = decltype(type);
        dataRepository::ViewWrapper<fieldType> *  to = elementSubRegion->template RegisterViewWrapper< fieldType >(m_to);
        dataRepository::ViewWrapper<fieldType> & from = dynamic_cast< dataRepository::ViewWrapper<fieldType> & >(*vwFrom);
        fieldType & fromField = from.reference();
        fieldType & toField = to->reference();
        toField = fromField;
      });
    });
  }

}

REGISTER_CATALOG_ENTRY( TaskBase, CopyField, std::string const &, ManagedGroup * const )

} /* namespace */
