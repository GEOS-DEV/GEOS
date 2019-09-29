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

#include "ElementRegionBase.hpp"

#include "common/TimingMacros.hpp"

//#include "constitutive/ConstitutiveManager.hpp"
//#include "finiteElement/FiniteElementDiscretizationManager.hpp"
//#include "finiteElement/basis/BasisBase.hpp"
//#include "finiteElement/quadrature/QuadratureBase.hpp"
//#include "managers/NumericalMethodsManager.hpp"
//#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;
//using namespace constitutive;


ElementRegionBase::ElementRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_numericalMethod()  //,
//    m_toNodesRelation(this->registerWrapper< array2d<integer>
// >(keys::nodeList).reference())
{
//  m_toNodesRelation.resize2(0,8);
//  this->registerWrapper<mapPair_array>(keys::constitutiveMap)->setSizedFromParent(1);

  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  this->RegisterGroup(viewKeyStruct::elementSubRegions);

  registerWrapper( viewKeyStruct::materialListString, &m_materialList, 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of materials present in this region");

}


ElementRegionBase::~ElementRegionBase()
{}


//REGISTER_CATALOG_ENTRY( ObjectManagerBase, ElementRegionBase, std::string const &, Group * const )

}
