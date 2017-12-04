/*
 * CellBlockSubRegion.hpp
 *
 *  Created on: May 11, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_

#include "CellBlock.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const constitutiveMap = "constitutiveMap";
string const constitutiveGrouping = "ConstitutiveGrouping";
}
}

class CellBlockSubRegion : public CellBlock
{
public:
  CellBlockSubRegion( string const & name, ManagedGroup * const parent );
  virtual ~CellBlockSubRegion();

  void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  void ReadXML_PostProcess() override;

  void InitializePreSubGroups( ManagedGroup * const ) override final;

  void InitializePostSubGroups( ManagedGroup * const ) override final;

  void CopyFromCellBlock( CellBlock const * source );

  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {
    auto const & constitutiveGrouping = this->getReference< map< string, localIndex_array > >(dataRepository::keys::constitutiveGrouping);

    for( auto & constitutiveGroup : constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  void MaterialPassThru( string const & matName,
                         string const & setName,
                         lSet & materialSet,
                         ManagedGroup * material );

  struct viewKeysStruct : public CellBlock::viewKeysStruct
  {} viewKeys;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_ */
