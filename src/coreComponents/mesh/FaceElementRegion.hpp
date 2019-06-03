/**
 * @file FaceElementRegion.hpp
 *
 */

#ifndef CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_
#define CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_

#include "ElementRegion.hpp"

namespace geosx
{

/**
 * @class FaceElementRegion
 *
 * The FaceElementRegion class contains the functionality to support the concept of a FaceElementRegion in the element
 * hierarchy. FaceElementRegion derives from ElementRegion and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class FaceElementRegion : public ElementRegion
{
public:
  /**
   * @brief constructor
   * @param name The name of the object in the data hierarchy.
   * @param parent Pointer to the parent group in the data hierarchy.
   */
  FaceElementRegion( string const & name, ManagedGroup * const parent );

  FaceElementRegion() = delete;
  virtual ~FaceElementRegion() override;

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string CatalogName()
  { return "FaceElementRegion"; }

  virtual const string getCatalogName() const override final
  { return FaceElementRegion::CatalogName(); }


  virtual void GenerateMesh( ManagedGroup const * ) override {}

  /**
   * @brief This function generates the face/fracture mesh
   * @param faceManager
   */
  void GenerateFractureMesh( FaceManager const * const faceManager );

  /**
   * @brief This function generates and adds entries to the face/fracture mesh
   * @param faceManager A pointer to the FaceManager object.
   * @param subRegionName The name of the FaceElementSubRegion to insert the new entries.
   * @param faceIndices The local indices of the new faces that define the face element.
   * @return The local index of the new FaceElement entry.
   */
  localIndex AddToFractureMesh( FaceManager const * const faceManager,
                                string const & subRegionName,
                                localIndex const faceIndices[2] );


  struct viewKeyStruct : public ElementRegion::viewKeyStruct
  {
    static constexpr auto fractureSetString = "fractureSet";
    static constexpr auto edgesTofractureConnectorsString = "edgesToFractureConnectors";
    static constexpr auto fractureConnectorsToEdgesString = "fractureConnectorsToEdges";
    static constexpr auto fractureConnectorsToFaceElementsString = "fractureElementConnectors";
    static constexpr auto faceElementsToCellsString = "fractureCellConnectors";
    static constexpr auto fractureCellConnectorIndicesString = "fractureCellConnectorIndices";
  };

  set< localIndex > m_recalculateConnectors;
  set< localIndex > m_newFractureElements;

private:
  string_array m_fractureSetNames;

  map< localIndex, localIndex > m_edgesToFractureConnectors;
  array1d<localIndex> m_fractureConnectorsToEdges;
  array1d< array1d<localIndex> > m_fractureConnectorsToFaceElements;
  array1d<localIndex > m_fractureCellConnectorIndices;
  FixedToManyElementRelation m_faceElementsToCells;

};

} /* namespace geosx */

#endif /* CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_ */
