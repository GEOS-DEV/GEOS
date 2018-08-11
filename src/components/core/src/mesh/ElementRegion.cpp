/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * ElementManagerT.cpp
 *
 *  Created on: Sep 14, 2010
 *      Author: settgast1
 */

#include "ElementRegion.hpp"
#include "CellBlockManager.hpp"
#include "CellBlockSubRegion.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/basis/BasisBase.hpp"
#include "finiteElement/quadrature/QuadratureBase.hpp"

#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;


ElementRegion::ElementRegion( string const & name, ManagedGroup * const parent ):
  ObjectManagerBase( name, parent )  //,
//    m_toNodesRelation(this->RegisterViewWrapper< array2d<integer>
// >(keys::nodeList).reference())
{
//  m_toNodesRelation.resize2(0,8);
//  this->RegisterViewWrapper<mapPair_array>(keys::constitutiveMap)->setSizedFromParent(1);
  this->RegisterGroup(keys::cellBlockSubRegions);
}


ElementRegion::~ElementRegion()
{}


void ElementRegion::FillDocumentationNode()
{
  ObjectManagerBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->getCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "an element region" );



  docNode->AllocateChildNode( keys::defaultMaterial,
                              keys::defaultMaterial,
                              -1,
                              "string",
                              "string",
                              "Default Material Name",
                              "Default Material Name",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::numericalMethod,
                              keys::numericalMethod,
                              -1,
                              "string",
                              "string",
                              "Default Material Name",
                              "Default Material Name",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );


//  docNode->AllocateChildNode( keys::constitutiveMap,
//                              keys::constitutiveMap,
//                              -1,
//                              "mapPair_array",
//                              "mapPair_array",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "NONE",
//                              "",
//                              1,
//                              0,
//                              0 );

  docNode->AllocateChildNode( keys::cellBlockSubRegionNames,
                              keys::cellBlockSubRegionNames,
                              -1,
                              "string_array",
                              "string_array",
                              "Number of Nodes Per Element",
                              "Number of Nodes Per Element",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );



}

void ElementRegion::ReadXML_PostProcess()
{
//  integer & numNodesPerElem = *(getData<integer>(keys::numNodesPerElement));
//  numNodesPerElem = 8;
}

/**
 *
 * @param problemManager
 * @param counts
 *
 * This function sets two mapping objects.
 *
 * The first is the constitutiveMap, which is a pair of 2D arrays. The size of
 * each array is ( numberOfElements x numberOfQuadraturePointsPerElement ). The
 * first array contains the material index.
 * The second array contains the index for a given material. This is a
 * preliminary implementation of this concept.
 * Ultimately, there will be an arbitrary number of material points per
 * quadrature point.
 *
 * The second mapping object is a so called material grouping array. It is a map
 * that contains the elements that contain
 * a given material. So the key is the material, and the value is an array of
 * element indices that contain that material.
 * Again, this is preliminary. This will be refined to contain elements that are
 * pure of a single material for all
 * quadrature points.
 *
 */
void ElementRegion::SetConstitutiveMap( ManagedGroup const * problemManager,
                                        map<string,localIndex> & counts )
{
//  map<string,integer> counts;
  ManagedGroup const * domain = problemManager->GetGroup(keys::domain);
  ConstitutiveManager const * constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

//  ConstitutiveManager::constitutiveMaps constitutiveMapPair =
// constitutiveManager->GetMaps( 1 );
  typename ManagedGroup::subGroupMap::LookupMapType const & constitutiveIndexLookup = constitutiveManager->GetSubGroups().keys();

  string defaultMaterial = this->getData<string>(keys::defaultMaterial);
  localIndex defaultMaterialIndex = constitutiveIndexLookup.at(defaultMaterial);


  auto const & numMethodName = this->getData<string>(keys::numericalMethod);
  FiniteElementManager const * numericalMethodManager = problemManager->GetGroup<FiniteElementManager>(keys::finiteElementManager);
  FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
  FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);
  auto const & quadratureName = feSpace->getData<string>(keys::quadrature);
  QuadratureBase const & quadrature = numericalMethodManager->GetGroup(keys::quadratureRules)->getReference<QuadratureBase>( quadratureName );


  ManagedGroup * cellBlockSubRegions = this->GetGroup(keys::cellBlockSubRegions);
  for( auto & cellSubBlock : cellBlockSubRegions->GetSubGroups() )
  {
    auto & cellToConstitutiveMap = cellSubBlock.second->getReference< std::pair< array2d< localIndex >, array2d< localIndex > > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString);
    auto & constitutiveGrouping = cellSubBlock.second->getReference< map< string, localIndex_array > >(CellBlockSubRegion::viewKeyStruct::constitutiveGroupingString);
//    constitutiveGrouping.resize( constitutiveMapPair.second.size() );
//    constitutiveGrouping["mix"];

//    localIndex counter = 0;
    for( localIndex k = 0 ; k < cellSubBlock.second->size() ; ++k )
    {
      constitutiveGrouping[defaultMaterial].push_back(k);
      for( localIndex q = 0 ; q < quadrature.size() ; ++q )
      {
        cellToConstitutiveMap.first(k,q)  = defaultMaterialIndex;
        cellToConstitutiveMap.second(k,q) = counts[defaultMaterial]++;
//        ++(counts[defaultMaterial]);
      }
    }
    // for( auto mat : constitutiveGrouping )
    // {
    //   for( auto a=0 ; a<mat.second.size() ; ++a )
    //   {
    //     std::cout<<cellSubBlock.second->getName()<<"
    // constitutiveGrouping["<<mat.first<<"]["<<a<<"] =
    // "<<mat.second[a]<<std::endl;
    //   }

    // }
  }
//  return counts;
}

void ElementRegion::InitializePreSubGroups( ManagedGroup * const problemManager )
{

  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>(keys::domain);
  ManagedGroup const * cellBlockManager = domain->GetGroup(keys::cellManager);

  ManagedGroup * cellBlockSubRegions = this->GetGroup(dataRepository::keys::cellBlockSubRegions);

  for( auto const & cellBlockName : this->getReference<string_array>(keys::cellBlockSubRegionNames) )
  {
    CellBlockSubRegion * cellBlock = cellBlockSubRegions->RegisterGroup<CellBlockSubRegion>(cellBlockName);
    cellBlock->FillDocumentationNode();
    cellBlock->RegisterDocumentationNodes();
  }

  auto const & numMethodName = this->getData<string>(keys::numericalMethod); 
  FiniteElementManager const * numericalMethodManager = problemManager->GetGroup<FiniteElementManager>(keys::finiteElementManager);
  FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
  FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

  auto const & basisName = feSpace->getData<string>(keys::basis);
  auto const & quadratureName = feSpace->getData<string>(keys::quadrature);
  BasisBase const & basis = numericalMethodManager->GetGroup(keys::basisFunctions)->getReference<BasisBase>( basisName );
  QuadratureBase const & quadrature = numericalMethodManager->GetGroup(keys::quadratureRules)->getReference<QuadratureBase>( quadratureName );

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  r1_array const & X = mesh->getNodeManager()->getReference<r1_array>(keys::referencePositionString);

  forCellBlocks([&]( CellBlockSubRegion * subRegion )
    {
      ManagedGroup const * cellBlocks = cellBlockManager->GetGroup(keys::cellBlocks);
      subRegion->CopyFromCellBlock( cellBlocks->GetGroup<CellBlock>( subRegion->getName() ) );

      feSpace->ApplySpaceToTargetCells(subRegion);


      feSpace->CalculateShapeFunctionGradients( X, subRegion);

//    feSpace
    });



}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ElementRegion, std::string const &, ManagedGroup * const )

}
