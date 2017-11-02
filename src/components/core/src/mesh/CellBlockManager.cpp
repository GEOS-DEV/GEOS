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

#include "CellBlockManager.hpp"

#include "FaceManager.hpp"
//#include "legacy/IO/BinStream.h"
#include <map>
#include <vector>
//#include "legacy/Constitutive/Material/MaterialFactory.h"
//#include "legacy/ArrayT/ArrayT.h"

namespace geosx
{
using namespace dataRepository;

CellBlockManager::CellBlockManager(  string const & name, ManagedGroup * const parent ):
ObjectManagerBase(name,parent)
{
  this->RegisterGroup<ManagedGroup>(keys::cellBlocks);
}

CellBlockManager::~CellBlockManager()
{
  // TODO Auto-generated destructor stub
}

void CellBlockManager::resize( integer_array const & numElements,
                             string_array const & regionNames,
                             string_array const & elementTypes )
{
  localIndex const numRegions = regionNames.size();
//  ManagedGroup * elementRegions = this->GetGroup(keys::cellBlocks);
  for( localIndex reg=0 ; reg<numRegions ; ++reg )
  {
    CellBlock * elemRegion = this->GetRegion( regionNames[reg] );
    elemRegion->resize(numElements[reg]);
  }
}


//CellBlock & CellBlockManager::CreateRegion( string const & regionName,
//                                             string const & elementType,
//                                             integer const & numElements )
//{
////  ElementRegion * elemRegion = elementRegions.RegisterGroup( regionNames );
////  elemRegion->resize(numElements);
//}

void CellBlockManager::CreateChild( string const & childKey, string const & childName )
{
}
//  ManagedGroup * elementRegions = this->GetGroup(keys::cellBlocks);
//  for (pugi::xml_node childNode=targetNode.first_child(); childNode; childNode=childNode.next_sibling())
//  {
//    if( childNode.name() == string("ElementRegion") )
//    {
//      std::string regionName = childNode.attribute("name").value();
//      std::cout<<regionName<<std::endl;
//
//      CellBlock * elemRegion = elementRegions->RegisterGroup<CellBlock>( regionName );
//      elemRegion->SetDocumentationNodes( nullptr );
//      elemRegion->RegisterDocumentationNodes();
//      elemRegion->ReadXML(childNode);
//    }
//  }
//}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlockManager, string const &, ManagedGroup * const )
}
