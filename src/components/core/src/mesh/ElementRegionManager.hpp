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
/**
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ZONEMANAGER_H
#define ZONEMANAGER_H

//#include "Common.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "CellBlock.hpp"
#include "CellBlockSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"

//#include "legacy/ArrayT/bufvector.h"
#include "ElementRegion.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const elementRegions = "elementRegions";
//string const elementRegionManager="ElementRegions";
}
}

/**
 * Class to manage the data stored at the element level.
 */
class ElementRegionManager : public ObjectManagerBase
{
public:

  template< typename VIEWTYPE >
  using ElementViewAccessor = array< array< VIEWTYPE * > > ;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "ZoneManager"; }

  virtual const string getCatalogName() const override final
  { return ElementRegionManager::CatalogName(); }



  ///@}

  ElementRegionManager( string const &, ManagedGroup * const parent );
  virtual ~ElementRegionManager() override;

  localIndex getNumberOfElements() const;

//  void Initialize(  ){}

  void InitializePreSubGroups( ManagedGroup * const ) override final;
  void InitializePostSubGroups( ManagedGroup * const ) override final;

  // virtual void CreateChild( string const & childKey, string const & childName ) override;
  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode ) override;


  using ManagedGroup::resize;

  void resize( integer_array const & numElements,
               string_array const & regionNames,
               string_array const & elementTypes );

//  CellBlock & CreateRegion( string const & regionName,
//                               string const & elementType,
//                               integer const & numElements );

  ElementRegion const * GetRegion( string const & regionName ) const
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(regionName);
  }
  ElementRegion * GetRegion( string const & regionName )
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(regionName);
  }

  ElementRegion const * GetRegion( localIndex const & index ) const
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(index);
  }
  ElementRegion * GetRegion( localIndex const & index )
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetGroup<ElementRegion>(index);
  }

  localIndex numRegions() const
  {
    return this->GetGroup(dataRepository::keys::elementRegions)->GetSubGroups().size();
  }

  localIndex numCellBlocks() const;


  template< typename LAMBDA >
  void forElementRegions( LAMBDA lambda )
  {
    ManagedGroup * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);
    elementRegions->forSubGroups<ElementRegion>( lambda );

  }

  template< typename LAMBDA >
  void forElementRegions( LAMBDA lambda ) const
  {
    ManagedGroup const * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);
    elementRegions->forSubGroups<ElementRegion>( lambda );
  }

  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda )
  {
    ManagedGroup * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);

    for( auto & region : elementRegions->GetSubGroups() )
    {
      ManagedGroup * cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
      for( auto & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
      {
        CellBlockSubRegion * cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);
        lambda( cellBlock );
      }
    }
  }

  template< typename LAMBDA >
  void forCellBlocks( LAMBDA lambda ) const
  {
    ManagedGroup const * elementRegions = this->GetGroup(dataRepository::keys::elementRegions);

    for( auto const & region : elementRegions->GetSubGroups() )
    {
      ManagedGroup const * cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
      for( auto const & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
      {
        CellBlockSubRegion const * cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);
        lambda( cellBlock );
      }
    }
  }

  template< typename VIEWTYPE >
  ElementViewAccessor<VIEWTYPE> ConstructViewAccessor( string const & name,
                                                       string const & neighborName );


  using ManagedGroup::PackSize;
  using ManagedGroup::Pack;
  using ObjectManagerBase::PackGlobalMapsSize;
  using ObjectManagerBase::PackGlobalMaps;
  using ObjectManagerBase::UnpackGlobalMaps;
  using ObjectManagerBase::PackUpDownMapsSize;
  using ObjectManagerBase::PackUpDownMaps;
  using ObjectManagerBase::UnpackUpDownMaps;



  int PackSize( array<string> const & wrapperNames,
                ElementViewAccessor<localIndex_array> const & packList ) const;

  int Pack( buffer_unit_type * & buffer,
            array<string> const & wrapperNames,
            ElementViewAccessor<localIndex_array> const & packList ) const;

  using ObjectManagerBase::Unpack;
  int Unpack( buffer_unit_type const * & buffer,
              ElementViewAccessor<localIndex_array> & packList );



  int PackGlobalMapsSize( ElementViewAccessor<localIndex_array> const & packList ) const;

  int PackGlobalMaps( buffer_unit_type * & buffer,
                              ElementViewAccessor<localIndex_array> const & packList ) const;


  int UnpackGlobalMaps( buffer_unit_type const * & buffer,
                                ElementViewAccessor<localIndex_array> & packList );

  int PackUpDownMapsSize( ElementViewAccessor<localIndex_array> const & packList ) const;

  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor<localIndex_array> const & packList ) const;


  int UnpackUpDownMaps( buffer_unit_type const * & buffer,
                        ElementViewAccessor<localIndex_array> const & packList );




private:
  template< bool DOPACK >
  int PackPrivate( buffer_unit_type * & buffer,
                   array<string> const & wrapperNames,
                   ElementViewAccessor<localIndex_array> const & viewAccessor ) const;

  template< bool DOPACK >
  int PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                             ElementViewAccessor<localIndex_array> const & viewAccessor ) const;

  template< bool DOPACK >
  int
  PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                         ElementViewAccessor<localIndex_array> const & packList ) const;

  ElementRegionManager( const ElementRegionManager& );
  ElementRegionManager& operator=( const ElementRegionManager&);
};



template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor<VIEWTYPE>
ElementRegionManager::
ConstructViewAccessor( string const & viewName,
                       string const & neighborName )
{
  ElementViewAccessor<VIEWTYPE> viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0 ; kReg<numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = GetRegion(kReg);
    viewAccessor[kReg].resize( elemRegion->numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(kSubReg);

      if( neighborName.empty() )
      {
        viewAccessor[kReg][kSubReg] = &(subRegion->getReference<VIEWTYPE>(viewName));
      }
      else
      {
        viewAccessor[kReg][kSubReg] = &(subRegion->GetGroup(ObjectManagerBase::groupKeyStruct::neighborDataString)->
                                      GetGroup(neighborName)->getReference<VIEWTYPE>(viewName));
      }
    }
  }

  return viewAccessor;

}





}
#endif /* ZONEMANAGER_H */
