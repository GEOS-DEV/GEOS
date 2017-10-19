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
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTREGION_H
#define ELEMENTREGION_H

#include "managers/ObjectManagerBase.hpp"
#include "legacy/ObjectManagers/EnergyT.h"
//#include "common/InterObjectRelation.hpp"
//#include "legacy/ArrayT/bufvector.h"
#include "FaceManager.hpp"



class StableTimeStep;

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const defaultMaterial = "material";
//string const numNodesPerElement = "numNodesPerElement";
//string const nodeList = "nodeList";
//string const constitutiveMap = "constitutiveMap";
string const numericalMethod = "numericalMethod";
string const cellBlockSubRegions = "cellBlockSubRegions";
string const cellBlockSubRegionNames = "cellBlocks";
}
}

class CellBlockSubRegion;



/**
 * Class to manage the data stored at the element level.
 */
class ElementRegion : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static string CatalogName()
  {
    return "ElementRegion";
  }

  string getCatalogName() const override final
  {
    return ElementRegion::CatalogName();
  }


  ///@}


  ElementRegion() = delete;

  ElementRegion( string const & name, ManagedGroup * const parent );


  ElementRegion(const ElementRegion& init);
//  ElementRegion( ElementRegion&& init);
  

  void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  virtual void ReadXML_PostProcess() override;

  void SetConstitutiveMap( dataRepository::ManagedGroup const * domain,
                                        map<string,localIndex> & counts );

  virtual ~ElementRegion();

//  Array2dT<integer> & m_toNodesRelation;


  virtual void InitializePreSubGroups( ManagedGroup * const group ) override final;


  CellBlockSubRegion const * GetSubRegion( string const & regionName ) const
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(regionName);
  }
  CellBlockSubRegion * GetSubRegion( string const & regionName )
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(regionName);
  }

  CellBlockSubRegion const * GetSubRegion( integer const & index ) const
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(index);
  }
  CellBlockSubRegion * GetSubRegion( integer const & index )
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetGroup<CellBlockSubRegion>(index);
  }

  integer numSubRegions() const
  {
    return this->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups().size();
  }
  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda )
  {
    ManagedGroup * cellBlockSubRegions = this->GetGroup(dataRepository::keys::cellBlockSubRegions);

    cellBlockSubRegions->forSubGroups<CellBlockSubRegion>( [&]( CellBlockSubRegion * subRegion ) -> void
    {
      lambda( subRegion );
    });

//    set<string> names;
//    for( auto & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
//    {
//      names.insert(iterCellBlocks.first);
//    }
//    for( auto & name : names )
//    {
//      lambda( cellBlockSubRegions->GetGroup<CellBlockSubRegion>(name) );
//    }

  }


  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda ) const
  {
    ManagedGroup const * cellBlockSubRegions = this->GetGroup(dataRepository::keys::cellBlockSubRegions);

    cellBlockSubRegions->forSubGroups<CellBlockSubRegion>( [&]( CellBlockSubRegion const * subRegion ) -> void
    {
      lambda( subRegion );
    });

//        set<string> names;
//    for( auto const & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
//    {
//      names.insert(iterCellBlocks.first);
//    }
//
//    for( auto const & name : names )
//    {
//      lambda( cellBlockSubRegions->GetGroup<CellBlockSubRegion>(name) );
//    }

  }


  string const & getNumericalMethod() const
  {
    return this->getData<string>(dataRepository::keys::numericalMethod);
  }


private:
  ElementRegion& operator=(const ElementRegion& rhs);

//  string & m_elementType;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////





}




#endif /* ELEMENTOBJECTT_H_ */
