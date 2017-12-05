/*
 * CellBlockSubRegion.cpp
 *
 *  Created on: May 11, 2017
 *      Author: rrsettgast
 */

#include "CellBlockSubRegion.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellBlockSubRegion::CellBlockSubRegion( string const & name, ManagedGroup * const parent ):
  CellBlock( name, parent )
{
  auto constitutiveGrouping = this->RegisterViewWrapper< map< string, localIndex_array > >(keys::constitutiveGrouping);
  constitutiveGrouping->setSizedFromParent(0);
}

CellBlockSubRegion::~CellBlockSubRegion()
{
  // TODO Auto-generated destructor stub
}



void CellBlockSubRegion::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->getCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "an element region" );

  CellBlock::FillDocumentationNode();

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

//  docNode->AllocateChildNode( keys::numNodesPerElement,
//                              keys::numNodesPerElement,
//                              -1,
//                              "integer",
//                              "integer",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              0,
//                              1,
//                              0 );



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



}

void CellBlockSubRegion::ReadXML_PostProcess()
{
//  integer & numNodesPerElem = numNodesPerElement();
//  numNodesPerElem = 8;

}

void CellBlockSubRegion::InitializePreSubGroups( ManagedGroup * const )
{
//  auto const & elementRegion = static_cast<ElementRegion const&>(
// *(this->getParent()) );
//  auto const & numMethod = elementRegion.getNumericalMethod();

}

void CellBlockSubRegion::InitializePostSubGroups( ManagedGroup * const )
{
  this->numNodesPerElement() = 8;
  this->numFacesPerElement() = 6;
}

void CellBlockSubRegion::CopyFromCellBlock( CellBlock const * source )
{
  this->resize(source->size());
  this->m_toNodesRelation = source->m_toNodesRelation;
}

void CellBlockSubRegion::MaterialPassThru( string const & matName,
                                           string const & setName,
                                           lSet & materialSet,
                                           ManagedGroup * material )
{}


} /* namespace geosx */
