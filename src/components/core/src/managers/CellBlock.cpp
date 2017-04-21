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

#include "CellBlock.hpp"

#include "Constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;



CellBlock::CellBlock( string const & name, ManagedGroup * const parent ):
    ObjectManagerBase( name, parent ),
    m_toNodesRelation(this->RegisterViewWrapper< Array2dT<real64> >(keys::nodeList).reference())
{
  m_toNodesRelation.resize2(0,8);
  this->RegisterViewWrapper<mapPair_array>(keys::constitutiveMap).setSizedFromParent(1);

}


CellBlock::~CellBlock()
{
}


void CellBlock::FillDocumentationNode( ManagedGroup * const group )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->getName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "an element region" );

  docNode->AllocateChildNode( keys::numNodesPerElement,
                              keys::numNodesPerElement,
                              -1,
                              "int32",
                              "int32",
                              "Number of Nodes Per Element",
                              "Number of Nodes Per Element",
                              "1",
                              "",
                              0,
                              1,
                              0 );
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
//                              0,
//                              1,
//                              0 );

//  docNode->AllocateChildNode( keys::numNodesPerElement,
//                              keys::numNodesPerElement,
//                              -1,
//                              "int32",
//                              "int32",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              1,
//                              0 );


}

void CellBlock::ReadXML_PostProcess()
{
  int32 & numNodesPerElem = *(getData<int32>(keys::numNodesPerElement));
  numNodesPerElem = 8;

}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlock, std::string const &, ManagedGroup * const )

}
