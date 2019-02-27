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

/*
 * ElementManagerT.cpp
 *
 *  Created on: Sep 14, 2010
 *      Author: settgast1
 */

#include "ElementRegion.hpp"

#include "CellBlockManager.hpp"
#include "CellElementSubRegion.hpp"
#include "FaceElementSubRegion.hpp"
#include "AggregateElementSubRegion.hpp"

#include "cxx-utilities/src/src/SparsityPattern.hpp"
#include "cxx-utilities/src/src/CRSMatrix.hpp"

#include "metis.h"
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


ElementRegion::ElementRegion( string const & name, ManagedGroup * const parent ):
  ObjectManagerBase( name, parent ),
  m_numericalMethod()  //,
//    m_toNodesRelation(this->RegisterViewWrapper< array2d<integer>
// >(keys::nodeList).reference())
{
//  m_toNodesRelation.resize2(0,8);
//  this->RegisterViewWrapper<mapPair_array>(keys::constitutiveMap)->setSizedFromParent(1);

  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  this->RegisterGroup(viewKeyStruct::elementSubRegions);

  RegisterViewWrapper( viewKeyStruct::materialListString, &m_materialList, 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of materials present in this region");

  RegisterViewWrapper( viewKeyStruct::sourceCellBlockNames, &m_cellBlockNames, false )->
    setInputFlag(InputFlags::OPTIONAL);

  RegisterViewWrapper( viewKeyStruct::fractureSetString, &m_fractureSetNames, false )->
    setInputFlag(InputFlags::OPTIONAL);

  RegisterViewWrapper( viewKeyStruct::coarseningRatioString, &m_coarseningRatio, false )->
    setInputFlag(InputFlags::OPTIONAL);

}


ElementRegion::~ElementRegion()
{}

void ElementRegion::PostProcessInput()
{
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
//void ElementRegion::SetConstitutiveMap( ManagedGroup const * problemManager,
//                                        map<string,localIndex> & counts )
//{
////  map<string,integer> counts;
//  ManagedGroup const * domain = problemManager->GetGroup(keys::domain);
//  ConstitutiveManager const * constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);
//
////  ConstitutiveManager::constitutiveMaps constitutiveMapPair =
//// constitutiveManager->GetMaps( 1 );
//  typename ManagedGroup::subGroupMap::LookupMapType const & constitutiveIndexLookup = constitutiveManager->GetSubGroups().keys();
//
//  string defaultMaterial = this->getData<string>(keys::defaultMaterial);
//  localIndex defaultMaterialIndex = constitutiveIndexLookup.at(defaultMaterial);
//
//
//  auto const & numMethodName = this->getData<string>(keys::numericalMethod);
//  NumericalMethodsManager const * numericalMethodManager = problemManager->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
//  FiniteElementSpaceManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
//  FiniteElementSpace const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementSpace>(numMethodName);
//  auto const & quadratureName = feDiscretization->getData<string>(keys::quadrature);
//  QuadratureBase const & quadrature = numericalMethodManager->GetGroup(keys::quadratureRules)->getReference<QuadratureBase>( quadratureName );
//
//
//  ManagedGroup * elementSubRegions = this->GetGroup(keys::cellBlockSubRegions);
//  for( auto & cellSubBlock : cellBlockSubRegions->GetSubGroups() )
//  {
//    auto & cellToConstitutiveMap = cellSubBlock.second->getReference< std::pair< array2d< localIndex >, array2d< localIndex > > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString);
//    auto & constitutiveGrouping = cellSubBlock.second->getReference< map< string, localIndex_array > >(CellBlockSubRegion::viewKeyStruct::constitutiveGroupingString);
////    constitutiveGrouping.resize( constitutiveMapPair.second.size() );
////    constitutiveGrouping["mix"];
//
////    localIndex counter = 0;
//    for( localIndex k = 0 ; k < cellSubBlock.second->size() ; ++k )
//    {
//      constitutiveGrouping[defaultMaterial].push_back(k);
//      for( localIndex q = 0 ; q < quadrature.size() ; ++q )
//      {
//        cellToConstitutiveMap.first(k,q)  = defaultMaterialIndex;
//        cellToConstitutiveMap.second(k,q) = counts[defaultMaterial]++;
////        ++(counts[defaultMaterial]);
//      }
//    }
//    // for( auto mat : constitutiveGrouping )
//    // {
//    //   for( auto a=0 ; a<mat.second.size() ; ++a )
//    //   {
//    //     std::cout<<cellSubBlock.second->getName()<<"
//    // constitutiveGrouping["<<mat.first<<"]["<<a<<"] =
//    // "<<mat.second[a]<<std::endl;
//    //   }
//
//    // }
//  }
////  return counts;
//}

//void ElementRegion::HangConstitutiveRelations( ManagedGroup const * problemManager )
//{
//  string const & numMethodName = this->getReference<string>(keys::numericalMethod);
//  NumericalMethodsManager const * numericalMethodManager = problemManager->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
//  int quadratureSize = 1;
//  ManagedGroup const * domain = problemManager->GetGroup(keys::domain);
//  ConstitutiveManager const * constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);
//  FiniteElementSpaceManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementDiscretizations);
//    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(numMethodName);
//  if( feDiscretization)
//  {
//    string const & quadratureName = feDiscretization->getReference<string>(keys::quadrature);
//    QuadratureBase const & quadrature = numericalMethodManager->GetGroup(keys::quadratureRules)->getReference<QuadratureBase>( quadratureName );
//    quadratureSize = quadrature.size() ;
//  }
//  forElementSubRegionsIndex( [&] ( localIndex const esr, CellBlockSubRegion * subRegion ) -> void
//      {
//      for( auto & materialName : m_materialList )
//      {
//        constitutiveManager->HangConstitutiveRelation( materialName, subRegion, quadratureSize );
//      }
//      });
//}

void ElementRegion::GenerateMesh( ManagedGroup const * const cellBlocks )
{
  ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);

  for( string const & cellBlockName : this->m_cellBlockNames )
  {
    CellElementSubRegion * subRegion = elementSubRegions->RegisterGroup<CellElementSubRegion>(cellBlockName);
    CellBlock const * source = cellBlocks->GetGroup<CellBlock>( subRegion->getName() );
    GEOS_ERROR_IF(source == nullptr, "Cell block named " + subRegion->getName() + " does not exist");
    subRegion->CopyFromCellBlock( source );
  }
}

void ElementRegion::GenerateAggregates( FaceManager const * const faceManager )
{
  std::cout << "Coarsening ratio is : " << m_coarseningRatio << std::endl;
  if(m_coarseningRatio == 0.)
  {
    return;
  }
  ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);
  AggregateElementSubRegion * const aggregateSubRegion =
    elementSubRegions->RegisterGroup<AggregateElementSubRegion>("coarse");

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();
  constexpr localIndex numElems = 2;

  // Counting the total number of cell and number of vertices  
  localIndex nbCellElements = 0;
  localIndex nbTotalVertices = 0;
  localIndex nbVertices = 0;
  this->forElementSubRegions( [&]( auto * const elementSubRegion ) -> void
    {
      nbCellElements += elementSubRegion->size();
      nbTotalVertices += elementSubRegion->size() * elementSubRegion->numNodesPerElement(0);
    });
  // Number of aggregate computation
  localIndex nbAggregates = integer_conversion< localIndex >( int(nbCellElements * m_coarseningRatio) );
  GEOS_LOG_RANK_0("Generating " +  std::to_string(nbAggregates) + " aggregates on region " + this->getName());

  // METIS variable declarations
  using idx_t = ::idx_t;
  idx_t options[METIS_NOPTIONS];                                    // Contains the METIS options
  METIS_SetDefaultOptions(options);                                 // ... That are set by default
  idx_t nnodes = integer_conversion< idx_t >( nbCellElements );     // Number of connectivity graph nodes
  idx_t nconst = 1;                                                 // Number of balancy constraints
  idx_t objval;                                                     // Total communication volume
  array1d< idx_t > parts(nnodes);                                // Map element index -> aggregate index
  //idx_t * parts = new idx_t[nnodes];                                // Map element index -> aggregate index
  idx_t nparts = integer_conversion< idx_t >( nbAggregates );       // Number of aggregates to be generated
  // METIS graph structure
  array1d< idx_t > xadjs(nnodes + 1);
  //idx_t * xadjs = new idx_t[nnodes+1];

  // Compute the connectivity graph
  LvArray::SparsityPattern< localIndex> graph(nbCellElements, nbCellElements);
  localIndex nbConnections = 0;
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    if( elemRegionList[kf][0] != -1 && elemRegionList[kf][1] != -1 )
    {
      localIndex const esr0 = elemSubRegionList[kf][0];
      localIndex const ei0  = elemList[kf][0];
      localIndex const esr1 = elemSubRegionList[kf][1];
      localIndex const ei1  = elemList[kf][1];
      graph.insertNonZero(ei0, ei1); //TODO : handle the case multiple region / element sub region
      graph.insertNonZero(ei1, ei0); //TODO : handle the case multiple region / element sub region
      nbConnections++;
    }
  }
  array1d< idx_t> adjncy(nbConnections*2);
  //idx_t * adjncy = new idx_t[integer_conversion< idx_t > (nbConnections*2)];
  
  // Fill the METIS graph structure
  xadjs[0] = integer_conversion< idx_t >(0);
  for( localIndex row = 0; row < nbCellElements; ++row )
  {
    localIndex numNonZeros = graph.numNonZeros( row );
    auto columns = graph.getColumns( row );
    xadjs[integer_conversion<idx_t>(row+1)] = xadjs[row] + integer_conversion< idx_t >( numNonZeros ) ;
    //xadjs[0] = 0 ;
    for(localIndex col = 0; col < numNonZeros; col++)
    {
      adjncy[xadjs[row]+col] = integer_conversion< idx_t >( columns[col] );
    }
  }
  
  // METIS partitionning
  METIS_PartGraphRecursive( &nnodes, &nconst, xadjs.data(), adjncy.data(), nullptr, nullptr, nullptr, &nparts,
                            nullptr, nullptr, options, &objval, parts.data());
  //aggregateSubRegion->CreateFromFineToCoarseMap(parts);
}

 void ElementRegion::GenerateFractureMesh( FaceManager const * const faceManager )
 {

   if( this->m_fractureSetNames.empty() )
   {
     return;
   }

   // key is edge index, value is faceElementIndex....this only works for a single fracture Region with a single subregion!!
   map< localIndex, set<localIndex> > fractureConnectorIndicesMap;

   array1d< localIndex > &
   fractureConnectorIndices = RegisterViewWrapper< array1d<localIndex > >( viewKeyStruct::fractureConnectorIndicesString )
     ->setRestartFlags( RestartFlags::NO_WRITE)
     ->setSizedFromParent(0)
     ->reference();

   array1d<array1d<localIndex> > &
   fractureConnectors = RegisterViewWrapper< array1d<array1d<localIndex> > >( viewKeyStruct::fractureElementConnectorString )
     ->setRestartFlags( RestartFlags::NO_WRITE)
     ->setSizedFromParent(0)
     ->reference();

   array1d< localIndex > &
   fractureCellConnectorIndices = RegisterViewWrapper< array1d<localIndex > >( viewKeyStruct::fractureCellConnectorIndicesString )
     ->setRestartFlags( RestartFlags::NO_WRITE)
     ->setSizedFromParent(0)
     ->reference();

   FixedToManyElementRelation &
   fractureCellConnectors = this->RegisterViewWrapper< FixedToManyElementRelation >( viewKeyStruct::fractureToCellConnectorString )
     ->setRestartFlags( RestartFlags::NO_WRITE)
     ->setSizedFromParent(0)
     ->reference();


   array2d<localIndex > const & faceToElementRegion = faceManager->elementRegionList();
   array2d<localIndex > const & faceToElementSubRegion = faceManager->elementSubRegionList();
   array2d<localIndex > const & faceToElementIndex = faceManager->elementList();

   ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);
   for( string const & setName : this->m_fractureSetNames )
   {
     FaceElementSubRegion * const subRegion = elementSubRegions->RegisterGroup<FaceElementSubRegion>(setName);
     set<localIndex> const & targetSet = faceManager->sets()->getReference<set<localIndex> >(setName);
     subRegion->resize( targetSet.size() );

     fractureCellConnectors.resize( targetSet.size(), 2 );

     FaceElementSubRegion::NodeMapType & nodeMap = subRegion->nodeList();
     FaceElementSubRegion::EdgeMapType & edgeMap = subRegion->edgeList();
     FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();

     OrderedVariableOneToManyRelation const & facesToNodesMap = faceManager->nodeList();
     OrderedVariableOneToManyRelation const & facesToEdgesMap = faceManager->edgeList();

     localIndex kfe = 0;
     std::cout << faceMap[0][0] << " " << faceMap[0][1] << std::endl;
     for( auto const faceIndex : targetSet )
     {
       faceMap[kfe][0] = faceIndex;
       faceMap[kfe][1] = faceIndex;

       arrayView1d<localIndex const> const & faceToNodesMap = facesToNodesMap[faceIndex];
       nodeMap[kfe].resize( faceToNodesMap.size() * 2 );
       for( localIndex a=0 ; a<faceToNodesMap.size() ; ++a )
       {
         const localIndex aa = a == 0 ? a : faceToNodesMap.size() - a;

         // TODO HACK need to generalize to something other than quads
         nodeMap[kfe][a] = faceToNodesMap[a];
         nodeMap[kfe][a+4] = faceToNodesMap[aa];
       }

       arrayView1d<localIndex const> const & faceToEdgesMap = facesToEdgesMap[faceIndex];
       edgeMap[kfe].resize( faceToEdgesMap.size() );
       for( localIndex a=0 ; a<faceToEdgesMap.size() ; ++a )
       {
         edgeMap[kfe][a] = faceToEdgesMap[a];
         fractureConnectorIndicesMap[ faceToEdgesMap[a] ].insert( kfe );
       }

       for( localIndex ke=0 ; ke<2 ; ++ke )
       {
         fractureCellConnectors.m_toElementRegion[kfe][ke] = faceToElementRegion[faceIndex][ke];
         fractureCellConnectors.m_toElementSubRegion[kfe][ke] = faceToElementSubRegion[faceIndex][ke];
         fractureCellConnectors.m_toElementIndex[kfe][ke] = faceToElementIndex[faceIndex][ke];
       }
       ++kfe;
     }
   }

   fractureConnectorIndices.resize( fractureConnectorIndicesMap.size() );
   fractureConnectors.resize( fractureConnectorIndicesMap.size() );
   localIndex connectorIndex=0;
   for( auto const & connector : fractureConnectorIndicesMap )
   {
     if( connector.second.size() > 1 )
     {
       fractureConnectorIndices[connectorIndex] = connector.first;
       fractureConnectors[connectorIndex].resize( connector.second.size() );
       localIndex fractureElementCounter = -1;
       for( auto const fractureElementIndex : connector.second )
       {
         fractureConnectors[connectorIndex][++fractureElementCounter] = fractureElementIndex;
       }
       ++connectorIndex;
     }
   }
   fractureConnectorIndices.resize(connectorIndex);
   fractureConnectors.resize(connectorIndex);


   forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion  * const subRegion )
   {
     FaceElementSubRegion::FaceMapType const & faceMap = subRegion->faceList();
     for( auto const & setIter : faceManager->sets()->wrappers() )
     {
       set<localIndex> const & faceSet = faceManager->sets()->getReference<set<localIndex> >( setIter.first );
       set<localIndex> & faceElementSet = subRegion->sets()->RegisterViewWrapper< set<localIndex> >( setIter.first )->reference();
       for( localIndex a=0 ; a<faceMap.size(0) ; ++a )
       {
         localIndex const faceIndex = faceMap[a][0];
         if( faceSet.count( faceIndex ) )
         {
           faceElementSet.insert( a );
         }
       }
     }
   });

 }


REGISTER_CATALOG_ENTRY( ObjectManagerBase, ElementRegion, std::string const &, ManagedGroup * const )

}
