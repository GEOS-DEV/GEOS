//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
#include "finiteElement/basis/BasisBase.hpp"
#include "finiteElement/quadrature/QuadratureBase.hpp"

#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;


ElementRegion::ElementRegion( string const & name, ManagedGroup * const parent ):
    ObjectManagerBase( name, parent )//,
//    m_toNodesRelation(this->RegisterViewWrapper< Array2dT<int32> >(keys::nodeList).reference())
{
//  m_toNodesRelation.resize2(0,8);
//  this->RegisterViewWrapper<mapPair_array>(keys::constitutiveMap)->setSizedFromParent(1);
  this->RegisterGroup(keys::cellBlockSubRegions);
}


ElementRegion::~ElementRegion()
{
}


void ElementRegion::FillDocumentationNode( ManagedGroup * const group )
{
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


  docNode->AllocateChildNode( keys::constitutiveMap,
                              keys::constitutiveMap,
                              -1,
                              "mapPair_array",
                              "mapPair_array",
                              "Number of Nodes Per Element",
                              "Number of Nodes Per Element",
                              "1",
                              "",
                              1,
                              0,
                              0 );

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
//  int32 & numNodesPerElem = *(getData<int32>(keys::numNodesPerElement));
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
 * each array is ( numberOfElements x numberOfQuadraturePointsPerElement ). The first array contains the material index.
 * The second array contains the index for a given material. This is a preliminary implementation of this concept.
 * Ultimately, there will be an arbitrary number of material points per quadrature point.
 *
 * The second mapping object is a so called material grouping array. It is a map that contains the elements that contain
 * a given material. So the key is the material, and the value is an array of element indices that contain that material.
 * Again, this is preliminary. This will be refined to contain elements that are pure of a single material for all
 * quadrature points.
 *
 */
void ElementRegion::SetConstitutiveMap( ManagedGroup const * problemManager,
                                                     map<string,localIndex> & counts )
{
//  map<string,int32> counts;
  ManagedGroup const * domain = problemManager->GetGroup(keys::domain);
  ConstitutiveManager const * constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

//  ConstitutiveManager::constitutiveMaps constitutiveMapPair = constitutiveManager->GetMaps( 1 );
  typename ManagedGroup::subGroupMap::LookupMapType const & constitutiveIndexLookup = constitutiveManager->GetSubGroups().keys();

  string defaultMaterial = this->getData<string>(keys::defaultMaterial);
  integer defaultMaterialIndex = constitutiveIndexLookup.at(defaultMaterial);


  auto const & numMethodName = this->getData<string>(keys::numericalMethod);
  FiniteElementManager const * numericalMethodManager = problemManager->GetGroup<FiniteElementManager>(keys::finiteElementManager);
  FiniteElementSpace const * feSpace = numericalMethodManager->GetGroup<FiniteElementSpace>(numMethodName);
  auto const & quadratureName = feSpace->getData<string>(keys::quadrature) ;
  QuadratureBase const & quadrature = numericalMethodManager->GetGroup(keys::quadratureRules)->getReference<QuadratureBase>( quadratureName );


  ManagedGroup * cellBlockSubRegions = this->GetGroup(keys::cellBlockSubRegions);
  for( auto & cellSubBlock : cellBlockSubRegions->GetSubGroups() )
  {
    auto & cellToConstitutiveMap = cellSubBlock.second->getReference< std::pair< Array2dT< integer >, Array2dT< integer > > >(keys::constitutiveMap);
    auto & constitutiveGrouping = cellSubBlock.second->getReference< map< string, integer_array > >(keys::constitutiveGrouping);
//    constitutiveGrouping.resize( constitutiveMapPair.second.size() );
//    constitutiveGrouping["mix"];
    
//    localIndex counter = 0;
    for( localIndex k = 0 ; k < cellSubBlock.second->size() ; ++k )
    {
      constitutiveGrouping[defaultMaterial].push_back(k);
      for( localIndex q = 0; q < quadrature.size(); ++q )
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
    //     std::cout<<cellSubBlock.second->getName()<<" constitutiveGrouping["<<mat.first<<"]["<<a<<"] = "<<mat.second[a]<<std::endl;
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
    cellBlock->FillDocumentationNode(nullptr);
    cellBlock->RegisterDocumentationNodes();
  }


  auto const & numMethodName = this->getData<string>(keys::numericalMethod);
  FiniteElementManager const * numericalMethodManager = problemManager->GetGroup<FiniteElementManager>(keys::finiteElementManager);
  FiniteElementSpace const * feSpace = numericalMethodManager->GetGroup<FiniteElementSpace>(numMethodName);
  auto const & basisName = feSpace->getData<string>(keys::basis) ;
  auto const & quadratureName = feSpace->getData<string>(keys::quadrature) ;
  BasisBase const & basis = numericalMethodManager->GetGroup(keys::basisFunctions)->getReference<BasisBase>( basisName );
  QuadratureBase const & quadrature = numericalMethodManager->GetGroup(keys::quadratureRules)->getReference<QuadratureBase>( quadratureName );

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  r1_array const & X = mesh->getNodeManager()->getReference<r1_array>(keys::ReferencePosition);

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
