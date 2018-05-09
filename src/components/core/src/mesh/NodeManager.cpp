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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * NodeManagerT.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "NodeManager.hpp"
//#include "managers/DomainPartition.hpp"
//#include "ObjectManagers/FaceManagerT.h"
//#include "ObjectManagers/EdgeManagerT.h"
//#include "ObjectManagers/ElementManagerT.h"
//#include "Utilities/Utilities.h"
//#include <fstream>
//#include "ElementRegionT.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace multidimensionalArray;

// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager( std::string const & name,
                          ManagedGroup * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition()
{
  RegisterViewWrapper(viewKeyStruct::referencePositionString, &m_referencePosition, false );
  this->RegisterViewWrapper< array<localIndex_array> >(viewKeyStruct::elementRegionListString);
  this->RegisterViewWrapper< array<localIndex_array> >(viewKeyStruct::elementSubRegionListString);
  this->RegisterViewWrapper< array<localIndex_array> >(viewKeyStruct::elementListString);

}



// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
/*
   NodeManagerT::NodeManagerT( const NodeManagerT& init ):
   ObjectDataStructureBaseT(init),
   DataLengths()(this->m_DataLengths),
   m_refposition(NULL),
   m_displacement(NULL),
   m_incrementalDisplacement(NULL),
   m_velocity(NULL),
   m_acceleration(NULL),
   m_force(NULL),
   m_mass(NULL),
   m_toElementsRelation(init.m_toElementsRelation),
   m_nodeToFaceMap(m_UnorderedVariableOneToManyMaps["nodeToFaceMap"]),
   m_nodeToEdgeMap(m_UnorderedVariableOneToManyMaps["nodeToEdgeMap"])
   {}
 */

// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::~NodeManager()
{}


//void NodeManager::Initialize()
//{
//  this->AddKeyedDataField<FieldInfo::referencePosition>();
//
//}



void NodeManager::FillDocumentationNode()
{
  ObjectManagerBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->getCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "a node manager" );

//  docNode->AllocateChildNode( keys::elementRegionMap,
//                              keys::elementRegionMap,
//                              -1,
//                              "integer_array",
//                              "integer_array",
//                              "map to element region",
//                              "map to element region",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );
//
//  docNode->AllocateChildNode( keys::elementSubRegionMap,
//                              keys::elementSubRegionMap,
//                              -1,
//                              "integer_array",
//                              "integer_array",
//                              "map to element sub regions",
//                              "map to element sub regions",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );
//
//  docNode->AllocateChildNode( keys::elementMap,
//                              keys::elementMap,
//                              -1,
//                              "localIndex_array",
//                              "localIndex_array",
//                              "map to element in a subregion",
//                              "map to element in a subregion",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );

  docNode->AllocateChildNode( keys::referencePositionString,
                              keys::referencePositionString,
                              -1,
                              "r1_array",
                              "r1_array",
                              "reference position of nodes",
                              "reference position of nodes",
                              "",
                              "",
                              1,
                              0,
                              0 );
}

void NodeManager::SetElementMaps( ElementRegionManager const * const elementRegionManager )
{
  array<localIndex_array> & elementRegionList = this->getReference< array<localIndex_array> >(viewKeyStruct::elementRegionListString);
  array<localIndex_array> & elementSubRegionList = this->getReference< array<localIndex_array> >(viewKeyStruct::elementSubRegionListString);
  array<localIndex_array> & elementList = this->getReference< array<localIndex_array> >(viewKeyStruct::elementListString);

  for( localIndex a=0 ; a<size() ; ++a )
  {
    elementRegionList[a].clear();
    elementSubRegionList[a].clear();
    elementList[a].clear();
  }

  for( typename dataRepository::indexType kReg=0 ; kReg<elementRegionManager->numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = elementRegionManager->GetRegion(kReg);

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);

      FixedOneToManyRelation const & elemsToNodes = subRegion->getReference<FixedOneToManyRelation>(subRegion->viewKeys.nodeList);

      for( localIndex ke=0 ; ke<subRegion->size() ; ++ke )
      {
        arrayView1d<localIndex const> const nodeList = elemsToNodes[ke];
        for( localIndex a=0 ; a<elemsToNodes.size(1) ; ++a )
        {
          localIndex nodeIndex = nodeList[a];
          elementRegionList[nodeIndex].push_back( kReg );
          elementSubRegionList[nodeIndex].push_back( kSubReg );
          elementList[nodeIndex].push_back( ke );
        }
      }
    }
  }
}


void NodeManager::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementSubRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementListString));
}


int NodeManager::PackUpDownMapsSize( localIndex_array const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

int NodeManager::PackUpDownMaps( buffer_unit_type * & buffer,
                             localIndex_array const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
int NodeManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                        localIndex_array const & packList ) const
{
  int packedSize = 0;


  return packedSize;
}


int NodeManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                               localIndex_array const & packList )
{
  int unPackedSize = 0;

  return unPackedSize;
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, std::string const &, ManagedGroup * const )

}
