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

#include "ElementRegion.hpp"

#include "CellBlockManager.hpp"
#include "CellElementSubRegion.hpp"
#include "FaceElementSubRegion.hpp"
#include "AggregateElementSubRegion.hpp"
#include "common/TimingMacros.hpp"
#include "cxx-utilities/src/src/ChaiVector.hpp"
#include "cxx-utilities/src/src/SparsityPattern.hpp"

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

void ElementRegion::GenerateAggregates( FaceManager const * const faceManager, NodeManager const * const nodeManager )
{
  GEOSX_MARK_FUNCTION;

  if(m_coarseningRatio <= 0.)
  {
    return;
  }
  ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);
  localIndex regionIndex = getIndexInParent();
  AggregateElementSubRegion * const aggregateSubRegion =
    elementSubRegions->RegisterGroup<AggregateElementSubRegion>("coarse");

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();

  constexpr localIndex numElems = 2;

  // Counting the total number of cell and number of vertices  
  localIndex nbCellElements = 0;
  this->forElementSubRegions( [&]( auto * const elementSubRegion ) -> void
    {
      nbCellElements += elementSubRegion->size();
    });
  // Number of aggregate computation
  localIndex nbAggregates = integer_conversion< localIndex >( int(nbCellElements * m_coarseningRatio) );
  GEOS_LOG_RANK_0("Generating " << nbAggregates  << " aggregates on region " << this->getName());

  // METIS variable declarations
  using idx_t = ::idx_t;
  idx_t options[METIS_NOPTIONS];                                    // Contains the METIS options
  METIS_SetDefaultOptions(options);                                 // ... That are set by default
  idx_t nnodes = integer_conversion< idx_t >( nbCellElements );     // Number of connectivity graph nodes
  idx_t nconst = 1;                                                 // Number of balancy constraints
  idx_t objval;                                                     // Total communication volume
  array1d< idx_t > parts(nnodes);                                   // Map element index -> aggregate index
  idx_t nparts = integer_conversion< idx_t >( nbAggregates );       // Number of aggregates to be generated
  

  // Compute the connectivity graph
  LvArray::SparsityPattern< idx_t, idx_t > graph( integer_conversion< idx_t >( nbCellElements ),
                                                  integer_conversion< idx_t >( nbCellElements ) );
  localIndex nbConnections = 0;
  array1d< localIndex > offsetSubRegions( this->GetSubRegions().size() );
  for( localIndex subRegionIndex = 1; subRegionIndex < offsetSubRegions.size(); subRegionIndex++ )
  {
    offsetSubRegions[subRegionIndex] = offsetSubRegions[subRegionIndex - 1] + this->GetSubRegion(subRegionIndex)->size();
  }
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    if( elemRegionList[kf][0] == regionIndex && elemRegionList[kf][1] == regionIndex && elemRegionList[kf][0] )
    {
      localIndex const esr0 = elemSubRegionList[kf][0];
      idx_t const ei0  = integer_conversion< idx_t >( elemList[kf][0] + offsetSubRegions[esr0] );
      localIndex const esr1 = elemSubRegionList[kf][1];
      idx_t const ei1  = integer_conversion< idx_t >( elemList[kf][1] + offsetSubRegions[esr1] );
      graph.insertNonZero(ei0, ei1);
      graph.insertNonZero(ei1, ei0);
      nbConnections++;
    }
  }

  // METIS partitionning
  idx_t * offsets = const_cast< idx_t* >( graph.getOffsets() );
  idx_t * columns = const_cast< idx_t* >( &graph.getColumns(0)[0] );
  METIS_PartGraphRecursive( &nnodes, &nconst, offsets, columns, nullptr, nullptr, nullptr,
                            &nparts, nullptr, nullptr, options, &objval, parts.data() );

  // Compute Aggregate barycenters
  array1d< R1Tensor > aggregateBarycenters( nparts );
  array1d< real64 > aggregateVolumes( nparts );
  array1d< real64 > normalizeVolumes( nbCellElements );
  
  // First, compute the volume of each aggregates
  this->forElementSubRegions( [&]( auto * const elementSubRegion ) -> void
  {
    localIndex const subRegionIndex = elementSubRegion->getIndexInParent();
    for(localIndex cellIndex = 0; cellIndex< elementSubRegion->size() ; cellIndex++)
    {
      if( elementSubRegion->GhostRank()[cellIndex] >= 0 )
        continue;
      aggregateVolumes[parts[cellIndex + offsetSubRegions[subRegionIndex]]] += elementSubRegion->getElementVolume()[cellIndex];
    }
  });

  // Second, compute the normalized volume of each fine elements
  this->forElementSubRegions( [&]( auto * const elementSubRegion ) -> void
  {
    localIndex const subRegionIndex = elementSubRegion->getIndexInParent();
    for(localIndex cellIndex = 0; cellIndex< elementSubRegion->size() ; cellIndex++)
    {
      if( elementSubRegion->GhostRank()[cellIndex] >= 0 )
        continue;
      normalizeVolumes[cellIndex + offsetSubRegions[subRegionIndex]] =
        elementSubRegion->getElementVolume()[cellIndex] / aggregateVolumes[parts[cellIndex + offsetSubRegions[subRegionIndex]]];
    }
  });

  // Third, normalize the centers
  this->forElementSubRegions( [&]( auto * const elementSubRegion ) -> void
  {
    localIndex const subRegionIndex = elementSubRegion->getIndexInParent();
    for(localIndex cellIndex = 0; cellIndex< elementSubRegion->size() ; cellIndex++)
    {
      if( elementSubRegion->GhostRank()[cellIndex] >= 0 )
        continue;
      R1Tensor center = elementSubRegion->getElementCenter()[cellIndex];
      center *= normalizeVolumes[cellIndex + offsetSubRegions[subRegionIndex]];
      aggregateBarycenters[parts[cellIndex + offsetSubRegions[subRegionIndex]]] += center;
    }
  });

  // Convert from metis to GEOSX types
  array1d< localIndex > partsGEOS( parts.size() );
  for( localIndex fineCellIndex = 0; fineCellIndex < partsGEOS.size(); fineCellIndex++ )
  {
    partsGEOS[fineCellIndex] = integer_conversion< localIndex >( parts[fineCellIndex] );
  }
  aggregateSubRegion->CreateFromFineToCoarseMap(nbAggregates, partsGEOS, aggregateBarycenters);
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, ElementRegion, std::string const &, ManagedGroup * const )

}
