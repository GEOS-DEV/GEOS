/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ElementRegionBase.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/solid/SolidBase.hpp"


namespace geosx
{
using namespace dataRepository;


ElementRegionBase::ElementRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_materialList(),
  m_meshBody()
{

  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->registerGroup( viewKeyStruct::elementSubRegions() );

  registerWrapper( viewKeyStruct::materialListString(), &m_materialList ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of materials present in this region" );

  registerWrapper( viewKeyStruct::meshBodyString(), &m_meshBody ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue("").
    setDescription( "Mesh body that contains this region" );

}


ElementRegionBase::~ElementRegionBase()
{}


void ElementRegionBase::postProcessInput()
{
  Group const & meshBody = this->getParent().getParent().getParent().getParent().getParent();
  string const & meshBodyName = meshBody.getName();
  Group const & meshBodies = meshBody.getParent();
  localIndex const numberOfMeshBodies = meshBodies.numSubGroups();

  if( numberOfMeshBodies == 1 )
  {
    if( m_meshBody=="" )
    {
      m_meshBody = meshBodyName;
    }
    GEOSX_ERROR_IF_EQ_MSG( m_meshBody,
                           meshBodyName,
                           "MeshBody specified does not match MeshBody in hierarchy.");
  }
  else
  {
    bool meshBodyFound = false;
    meshBodies.forSubGroups( [&] ( Group const & meshBody )
    {
      if( meshBody.getName()==m_meshBody )
      {
        meshBodyFound = true;
      }
    });
    GEOSX_ERROR_IF( !meshBodyFound, "MeshBody was not found");
  }


}


}
