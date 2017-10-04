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
/**
 * @file NodeManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */


#ifndef NODEMANAGERT_H_
#define NODEMANAGERT_H_

#include "ObjectManagerBase.hpp"
#include <string.h>
#include "CellBlockManager.hpp"


// *********************************************************************************************************************
// *********************************************************************************************************************
class SiloFile;

namespace geosx
{

class CellBlock;
class FaceManager;
class EdgeManager;

namespace dataRepository
{
namespace keys
{
std::string const nodeManager    = "NodeManager";
std::string const elementRegionMap("elementRegionMap");
std::string const elementSubRegionMap("elementSubRegionMap");
std::string const elementMap("elementMap");

}
}

using namespace dataRepository;
/**
 * @author Randolph Settgast
 *
 * The NodeManagerT class manages the node data using the ObjectDataStructureBaseT as a data manager.
 * This means that each field is stored in an array where each array entry corresponds to a node.
 */
class NodeManager : public ObjectManagerBase
{
public:


  /// default constructor
  NodeManager( std::string const & name,
               ManagedGroup * const parent );



  /// default destructor
  ~NodeManager();

  static string CatalogName() { return dataRepository::keys::nodeManager; }
  string getCatalogName() const      { return NodeManager::CatalogName(); }


  void FillDocumentationNode( ManagedGroup * const group ) override final;


//  void Initialize();




public:



  /** @name Maps
   * The Maps
   */
  ///@{


//  UnorderedVariableOneToManyRelation&  m_nodeToFaceMap;
//  UnorderedVariableOneToManyRelation&  m_nodeToEdgeMap;

//  UnorderedVariableOneToManyRelation& m_toCrackSurfacesRelation;

  ///@}


  struct viewKeysStruct
  {
    dataRepository::ViewKey referencePosition = { "ReferencePosition" };
    dataRepository::ViewKey totalDisplacement = { "TotalDisplacement" };
    dataRepository::ViewKey nodeList           = { "nodeList" };
    dataRepository::ViewKey numFacesPerElement = { "numFacesPerElement" };
    dataRepository::ViewKey faceList           = { "faceList" };

  }viewKeys;


  struct groupKeysStruct
  {
  }groupKeys;

  view_rtype_const<r1_array> referencePosition() const { return this->getData<r1_array>(viewKeys.referencePosition); }
  view_rtype<r1_array>       referencePosition()       { return this->getData<r1_array>(viewKeys.referencePosition); }
  view_rtype_const<r1_array> totalDisplacement() const { return this->getData<r1_array>(viewKeys.totalDisplacement); }
  view_rtype<r1_array>       totalDisplacement()       { return this->getData<r1_array>(viewKeys.totalDisplacement); }
protected:

private:
  /// copy constructor
  NodeManager() = delete;
  NodeManager( const NodeManager& init ) = delete;
  NodeManager& operator=( const NodeManager&) = delete;


};
// *********************************************************************************************************************
// *********************************************************************************************************************


// *********************************************************************************************************************


}


#endif /* NODEMANAGERT_H_ */
