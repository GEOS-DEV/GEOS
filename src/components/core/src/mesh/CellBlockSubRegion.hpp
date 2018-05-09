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

class CellBlockSubRegion : public CellBlock
{
public:
  CellBlockSubRegion( string const & name, ManagedGroup * const parent );
  virtual ~CellBlockSubRegion() override;

  void FillDocumentationNode() override final;

  void ReadXML_PostProcess() override;

  void InitializePreSubGroups( ManagedGroup * const ) override final;

  void InitializePostSubGroups( ManagedGroup * const ) override final;

  void CopyFromCellBlock( CellBlock const * source );

  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {

    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  void MaterialPassThru( string const & matName,
                         string const & setName,
                         lSet & materialSet,
                         ManagedGroup * material );


  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual int PackUpDownMapsSize( localIndex_array const & packList ) const override;

  virtual int PackUpDownMaps( buffer_unit_type * & buffer,
                              localIndex_array const & packList ) const override;

  virtual int UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                localIndex_array const & packList ) override;


  struct viewKeyStruct : public CellBlock::viewKeyStruct
  {
    static constexpr auto constitutiveGroupingString = "ConstitutiveGrouping";
    static constexpr auto constitutiveMapString = "ConstitutiveMap";
    static constexpr auto dNdXString = "dNdX";

    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString };
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString };
    dataRepository::ViewKey dNdX                  = { dNdXString };

  } viewKeys;


  map< string, localIndex_array > m_constitutiveGrouping;
  std::pair< Array2dT< localIndex >, Array2dT< localIndex > > m_constitutiveMapView;
  array< Array2dT<R1Tensor> > m_dNdX;

private:
  template< bool DOPACK >
  int PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                             localIndex_array const & packList ) const;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_CELLBLOCKSUBREGION_HPP_ */
