/*
 * FaceElementRegion.hpp
 *
 *  Created on: May 15, 2019
 *      Author: settgast
 */

#ifndef CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_
#define CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_

#include "ElementRegion.hpp"

namespace geosx
{

class FaceElementRegion : public ElementRegion
{
public:
  FaceElementRegion( string const & name, ManagedGroup * const parent );
  FaceElementRegion() = delete;
  virtual ~FaceElementRegion() override;

  static const string CatalogName()
  { return "FaceElementRegion"; }

  virtual const string getCatalogName() const override final
  { return FaceElementRegion::CatalogName(); }

  virtual void GenerateMesh( ManagedGroup const * ) override {}

  void GenerateFractureMesh( FaceManager const * const faceManager );

  localIndex AddToFractureMesh( FaceManager const * const faceManager, localIndex const faceIndices[2] );


  struct viewKeyStruct : public ElementRegion::viewKeyStruct
  {
    static constexpr auto fractureSetString = "fractureSet";
    static constexpr auto edgesTofractureConnectorsMapString = "edgesToFractureConnectors";
    static constexpr auto fractureConnectorToEdgeMapString = "fractureConnectorsToEdges";
    static constexpr auto fractureElementConnectorString = "fractureElementConnectors";
    static constexpr auto fractureToCellConnectorString = "fractureCellConnectors";
    static constexpr auto fractureCellConnectorIndicesString = "fractureCellConnectorIndices";
  };

private:
  string_array m_fractureSetNames;

  map< localIndex, localIndex > m_edgesToFractureConnectors;
  array1d<localIndex> m_fractureConnectorsToEdges;
  array1d< array1d<localIndex> > m_fractureElementConnectors;
  array1d<localIndex > m_fractureCellConnectorIndices;
  FixedToManyElementRelation m_fractureToCellConnectors;
};

} /* namespace geosx */

#endif /* CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_ */
