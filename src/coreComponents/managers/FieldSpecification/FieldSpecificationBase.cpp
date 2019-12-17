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

#include "FieldSpecificationBase.hpp"

#include "mpiCommunications/MpiWrapper.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/DomainPartition.hpp"


namespace geosx
{
using namespace dataRepository;

FieldSpecificationBase::FieldSpecificationBase( string const & name, Group * parent ):
  Group( name, parent ),
  m_normalizeBySetSize(false)
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  registerWrapper( viewKeyStruct::setNamesString, &m_setNames, 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("Name of sets that boundary condition is applied to.");

  registerWrapper( viewKeyStruct::objectPathString, &m_objectPath, 0 )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Path to the target field");

  registerWrapper( viewKeyStruct::fieldNameString, &m_fieldName, 0 )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of field that boundary condition is applied to.");

  registerWrapper( viewKeyStruct::componentString, &m_component, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Component of field (if tensor) to apply boundary condition to");

  registerWrapper( viewKeyStruct::directionString, &m_direction, 0 )->
//      setApplyDefaultValue(0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Direction to apply boundary condition to");

  registerWrapper( viewKeyStruct::functionNameString, &m_functionName, 0 )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of function that specifies variation of the BC");

  registerWrapper( viewKeyStruct::bcApplicationTableNameString, &m_bcApplicationFunctionName, 0 )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of table that specifies the on/off application of the bc.");

  registerWrapper( viewKeyStruct::scaleString, &m_scale, 0 )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Scale factor for value of BC.");


  registerWrapper( viewKeyStruct::initialConditionString, &m_initialCondition, 0 )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("BC is applied as an initial condition.");

  registerWrapper( viewKeyStruct::beginTimeString, &m_beginTime, 0 )->
    setApplyDefaultValue(-1.0e99)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("time at which BC will start being applied.");

  registerWrapper( viewKeyStruct::endTimeString, &m_endTime, 0 )->
    setApplyDefaultValue(1.0e99)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("time at which bc will stop being applied");

  
  registerWrapper( viewKeyStruct::setSizeScalingFactorString, &m_setSizeScalingFactor, 0 )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("size of the set on which the boundary condition is applied");
  
}


FieldSpecificationBase::~FieldSpecificationBase()
{}

FieldSpecificationBase::CatalogInterface::CatalogType&
FieldSpecificationBase::GetCatalog()
{
  static FieldSpecificationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void FieldSpecificationBase::InitializePreSubGroups( Group * const rootGroup )
{
  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  integer localSetSize = 0;

  // compute the local size of the set
  real64 dummyTime = 0.0;
  fsManager.Apply( dummyTime, domain, m_objectPath, m_fieldName,
                    [&]( FieldSpecificationBase const * const GEOSX_UNUSED_ARG( fs ),
                    string const &,
                    set<localIndex> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {
    integer_array& ghostRank = subRegion->getReference<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );    
    
    for( auto a : lset )
    {
      if (ghostRank[a] < 0)
      {
        ++localSetSize; 
      }
    }
 
  });

  // compute the global set size
  integer globalSetSize = 0;
  MpiWrapper::allReduce( &localSetSize, &globalSetSize, 1, MPI_SUM, MPI_COMM_GEOSX );

  // set the scaling factor
  m_setSizeScalingFactor = globalSetSize >= 1 ? 1.0 / globalSetSize : 1;
}

  
void FieldSpecificationBase::PostProcessInput()
{}


REGISTER_CATALOG_ENTRY( FieldSpecificationBase, FieldSpecificationBase, string const &, Group * const )

}
