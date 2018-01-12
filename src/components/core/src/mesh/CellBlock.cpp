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
 * ElementManagerT.cpp
 *
 *  Created on: Sep 14, 2010
 *      Author: settgast1
 */

#include "CellBlock.hpp"

#include "NodeManager.hpp"

namespace geosx
{
using namespace dataRepository;
//using namespace constitutive;


CellBlock::CellBlock( string const & name, ManagedGroup * const parent ):
  ObjectManagerBase( name, parent ),
  viewKeys(),
  groupKeys(),
  m_toNodesRelation(this->RegisterViewWrapper< Array2dT<localIndex> >(viewKeys.nodeList.Key())->reference()),
  m_toFacesRelation(this->RegisterViewWrapper< Array2dT<localIndex> >(viewKeys.faceList.Key())->reference())

{

  m_toNodesRelation.resize(0,8);
  m_toFacesRelation.resize(0,6);
//  this->RegisterViewWrapper<mapPair_array>(keys::constitutiveMap).setSizedFromParent(1);

}


CellBlock::~CellBlock()
{}


void CellBlock::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  ObjectManagerBase::FillDocumentationNode();

  docNode->setName( this->getCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "an element region" );

  docNode->AllocateChildNode( viewKeys.numNodesPerElement.Key(),
                              viewKeys.numNodesPerElement.Key(),
                              -1,
                              "integer",
                              "integer",
                              "Number of Nodes Per Element",
                              "Number of Nodes Per Element",
                              "1",
                              "",
                              0,
                              1,
                              0 );

//  docNode->AllocateChildNode( viewKeys.nodeList.Key(),
//                              viewKeys.nodeList.Key(),
//                              -1,
//                              "integer_array",
//                              "integer_array",
//                              "nodelist",
//                              "nodelist",
//                              "8",
//                              "",
//                              0,
//                              1,
//                              0 );

  docNode->AllocateChildNode( viewKeys.numFacesPerElement.Key(),
                              viewKeys.numFacesPerElement.Key(),
                              -1,
                              "integer",
                              "integer",
                              "Number of Faces Per Element",
                              "Number of Faces Per Element",
                              "6",
                              "",
                              0,
                              1,
                              0 );

//  docNode->AllocateChildNode( keys::defaultMaterial,
//                              keys::defaultMaterial,
//                              -1,
//                              "string",
//                              "string",
//                              "Default Material Name",
//                              "Default Material Name",
//                              "REQUIRED",
//                              "",
//                              0,
//                              1,
//                              0 );
//
//
//  docNode->AllocateChildNode( keys::constitutiveMap,
//                              keys::constitutiveMap,
//                              -1,
//                              "mapPair_array",
//                              "mapPair_array",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              1,
//                              0,
//                              0 );

//  docNode->AllocateChildNode( keys::numNodesPerElement,
//                              keys::numNodesPerElement,
//                              -1,
//                              "integer",
//                              "integer",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              1,
//                              0 );


}

void CellBlock::ReadXML_PostProcess()
{
//  integer & numNodesPerElem = this->numNodesPerElement();
//  numNodesPerElem = 8;
  this->numNodesPerElement() = 8;
  this->numFacesPerElement() = 6;

}

//map<string,integer> CellBlock::SetConstitutiveMap( ManagedGroup const * domain
// )
//{
//  map<string,integer> counts;
//  view_rtype<mapPair_array> cellToConstitutiveMap =
// this->getData<mapPair_array>(keys::constitutiveMap);
//  ConstitutiveManager const * constitutiveManager =
// domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);
//
//  ConstitutiveManager::constitutiveMaps constitutiveMapPair =
// constitutiveManager->GetMaps( 1 );
//
//  string defaultMaterial = this->getData<string>(keys::defaultMaterial);
//  integer defaultMaterialIndex =
// constitutiveMapPair.second.at(defaultMaterial);
//
//
//  localIndex counter = 0;
//  for( localIndex k=0 ; k<this->size() ; ++k )
//  {
//    cellToConstitutiveMap[k] = std::make_pair( defaultMaterialIndex, counter++
// );
//    ++(counts.at(defaultMaterial));
//  }
//  return counts;
//}


void CellBlock::GetFaceNodes( const localIndex elementIndex,
                              const localIndex localFaceIndex,
                              localIndex_array& nodeIndicies) const
{
  // get nodelist for this element
  const localIndex* const elemToNodeMap = m_toNodesRelation[elementIndex];

  // resize the nodeIndicies based on element type (this is wrong for some types
  // of elements)
  nodeIndicies.resize(4);

//  if (!m_elementGeometryID.compare(0, 4, "C3D8"))
  {
    if (localFaceIndex == 0)
    {
      nodeIndicies[0] = elemToNodeMap[0];
      nodeIndicies[1] = elemToNodeMap[1];
      nodeIndicies[2] = elemToNodeMap[5];
      nodeIndicies[3] = elemToNodeMap[4];
    }
    else if (localFaceIndex == 1)
    {
      nodeIndicies[0] = elemToNodeMap[0];
      nodeIndicies[1] = elemToNodeMap[2];
      nodeIndicies[2] = elemToNodeMap[3];
      nodeIndicies[3] = elemToNodeMap[1];
    }
    else if (localFaceIndex == 2)
    {
      nodeIndicies[0] = elemToNodeMap[0];
      nodeIndicies[1] = elemToNodeMap[4];
      nodeIndicies[2] = elemToNodeMap[6];
      nodeIndicies[3] = elemToNodeMap[2];
    }
    else if (localFaceIndex == 3)
    {
      nodeIndicies[0] = elemToNodeMap[1];
      nodeIndicies[1] = elemToNodeMap[3];
      nodeIndicies[2] = elemToNodeMap[7];
      nodeIndicies[3] = elemToNodeMap[5];
    }
    else if (localFaceIndex == 4)
    {
      nodeIndicies[0] = elemToNodeMap[3];
      nodeIndicies[1] = elemToNodeMap[2];
      nodeIndicies[2] = elemToNodeMap[6];
      nodeIndicies[3] = elemToNodeMap[7];
    }
    else if (localFaceIndex == 5)
    {
      nodeIndicies[0] = elemToNodeMap[4];
      nodeIndicies[1] = elemToNodeMap[5];
      nodeIndicies[2] = elemToNodeMap[7];
      nodeIndicies[3] = elemToNodeMap[6];
    }

  }
//  else if (!m_elementGeometryID.compare(0, 4, "C3D6"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//      nodeIndicies[2] = elemToNodeMap[5];
//      nodeIndicies[3] = elemToNodeMap[4];
//    }
//    else if (localFaceIndex == 1)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[2];
//      nodeIndicies[2] = elemToNodeMap[3];
//      nodeIndicies[3] = elemToNodeMap[1];
//    }
//    else if (localFaceIndex == 2)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[2];
//      nodeIndicies[2] = elemToNodeMap[4];
//      nodeIndicies[3] = std::numeric_limits<localIndex>::max();
//    }
//    else if (localFaceIndex == 3)
//    {
//      nodeIndicies[0] = elemToNodeMap[1];
//      nodeIndicies[1] = elemToNodeMap[3];
//      nodeIndicies[2] = elemToNodeMap[5];
//      nodeIndicies[3] = std::numeric_limits<localIndex>::max();
//    }
//    else if (localFaceIndex == 4)
//    {
//      nodeIndicies[0] = elemToNodeMap[2];
//      nodeIndicies[1] = elemToNodeMap[3];
//      nodeIndicies[2] = elemToNodeMap[5];
//      nodeIndicies[3] = elemToNodeMap[4];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 4, "C3D4"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[2];
//      nodeIndicies[2] = elemToNodeMap[1];
//    }
//    else if (localFaceIndex == 1)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//      nodeIndicies[2] = elemToNodeMap[3];
//    }
//    else if (localFaceIndex == 2)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[3];
//      nodeIndicies[2] = elemToNodeMap[2];
//    }
//    else if (localFaceIndex == 3)
//    {
//      nodeIndicies[0] = elemToNodeMap[1];
//      nodeIndicies[1] = elemToNodeMap[2];
//      nodeIndicies[2] = elemToNodeMap[3];
//    }
//  }
//
//  else if ( !m_elementGeometryID.compare(0,4,"CPE2") )
//  {
//    if( localFaceIndex == 0 )
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//    }
//  }
//
//  else if ( !m_elementGeometryID.compare(0,4,"CPE3") )
//  {
//    if( localFaceIndex == 0 )
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//    }
//    else if( localFaceIndex == 1 )
//    {
//      nodeIndicies[0] = elemToNodeMap[1];
//      nodeIndicies[1] = elemToNodeMap[2];
//    }
//    else if( localFaceIndex == 2 )
//    {
//      nodeIndicies[0] = elemToNodeMap[2];
//      nodeIndicies[1] = elemToNodeMap[0];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 4, "CPE4"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//    }
//    else if (localFaceIndex == 1)
//    {
//      nodeIndicies[0] = elemToNodeMap[1];
//      nodeIndicies[1] = elemToNodeMap[3];
//    }
//    else if (localFaceIndex == 2)
//    {
//      nodeIndicies[0] = elemToNodeMap[3];
//      nodeIndicies[1] = elemToNodeMap[2];
//    }
//    else if (localFaceIndex == 3)
//    {
//      nodeIndicies[0] = elemToNodeMap[2];
//      nodeIndicies[1] = elemToNodeMap[0];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 4, "STRI"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//    }
//    else if (localFaceIndex == 1)
//    {
//      nodeIndicies[0] = elemToNodeMap[1];
//      nodeIndicies[1] = elemToNodeMap[2];
//    }
//    else if (localFaceIndex == 2)
//    {
//      nodeIndicies[0] = elemToNodeMap[2];
//      nodeIndicies[1] = elemToNodeMap[0];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 3, "S4R"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//      nodeIndicies[2] = elemToNodeMap[2];
//      nodeIndicies[3] = elemToNodeMap[3];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 4, "TRSH"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = elemToNodeMap[0];
//      nodeIndicies[1] = elemToNodeMap[1];
//      nodeIndicies[2] = elemToNodeMap[2];
//    }
//  }
//
//  else
//  {
//    GEOS_ERROR("Error.  Don't know what kind of element this is and cannot
// build faces.");
//  }

}

R1Tensor CellBlock::GetElementCenter(localIndex k, const NodeManager& nodeManager, const bool useReferencePos) const
{

  view_rtype_const<r1_array> X = nodeManager.referencePosition();
//  view_rtype_const<r1_array> u = nodeManager.totalDisplacement();
  const localIndex* const nodelist = m_toNodesRelation[k];
  R1Tensor elementCenter(0.0);
  for ( localIndex a = 0 ; a < numNodesPerElement() ; ++a)
  {
    const localIndex b = nodelist[a];
    elementCenter += X[b];
    if(!useReferencePos)
      elementCenter += X[b];
  }
  elementCenter /= numNodesPerElement();

  return elementCenter;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlock, std::string const &, ManagedGroup * const )

}
