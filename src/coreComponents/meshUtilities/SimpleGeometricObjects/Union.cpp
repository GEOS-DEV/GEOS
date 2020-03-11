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

/*
 * @file Union.cpp
 * @brief Union of an arbitrary number of other geometric objects.
 * @param Objects Names of objects to unite
 */

#include "Union.hpp"

namespace geosx
{
using namespace dataRepository;

Union::Union( const std::string& name, Group * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_subObjects()
{
  registerWrapper( viewKeyStruct::subObjectsString, &m_subObjects, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of geometric objects to unite");
}

Union::~Union()
{}


bool Union::IsCoordInObject( const R1Tensor& coord ) const
{
  bool rval = false;
  dataRepository::Group const * geometries = this->getParent();

  for (auto subObjectName : m_subObjects)
  {
    SimpleGeometricObjectBase const * const subObject = geometries->GetGroup<SimpleGeometricObjectBase>(subObjectName);
    
    GEOSX_ERROR_IF(subObject==nullptr, "Object not found: " << subObjectName);

    if (subObject->IsCoordInObject(coord))
    {
      rval = true;
    }
  }
 
  return rval;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Union, std::string const &, Group * const )

} /* namespace geosx */
