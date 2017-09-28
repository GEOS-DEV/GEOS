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
 * NodeManagerT.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "NodeManager.hpp"
#include "DomainPartition.hpp"
//#include "ObjectManagers/FaceManagerT.h"
//#include "ObjectManagers/EdgeManagerT.h"
//#include "ObjectManagers/ElementManagerT.h"
//#include "Utilities/Utilities.h"
//#include <fstream>
//#include "ElementRegionT.hpp"
namespace geosx
{
using namespace dataRepository;

// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager( std::string const & name,
                          ManagedGroup * const parent ):
ObjectManagerBase( name, parent )
{

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
{
}


//void NodeManager::Initialize()
//{
//  this->AddKeyedDataField<FieldInfo::referencePosition>();
//
//}



void NodeManager::FillDocumentationNode( ManagedGroup * const group )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->getCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "a node manager" );

  docNode->AllocateChildNode( keys::elementRegionMap,
                              keys::elementRegionMap,
                              -1,
                              "int32_array",
                              "int32_array",
                              "map to element region",
                              "map to element region",
                              "",
                              "",
                              0,
                              0,
                              0 );

  docNode->AllocateChildNode( keys::elementSubRegionMap,
                              keys::elementSubRegionMap,
                              -1,
                              "int32_array",
                              "int32_array",
                              "map to element sub regions",
                              "map to element sub regions",
                              "",
                              "",
                              0,
                              0,
                              0 );

  docNode->AllocateChildNode( keys::elementMap,
                              keys::elementMap,
                              -1,
                              "int64_array",
                              "int64_array",
                              "map to element in a subregion",
                              "map to element in a subregion",
                              "",
                              "",
                              0,
                              0,
                              0 );

}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, std::string const &, ManagedGroup * const )

}

