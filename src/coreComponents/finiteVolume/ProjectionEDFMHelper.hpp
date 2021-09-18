#ifndef GEOSX_FINITEVOLUME_PROJECTIONEDFMHELPER_HPP_
#define GEOSX_FINITEVOLUME_PROJECTIONEDFMHELPER_HPP_

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{

class MeshLevel;
class CellElementStencilTPFA;
class EmbeddedSurfaceToCellStencil;

/**
 * @struct CellDescriptor
 * @brief A structure containing a single cell (element) identifier triplet.
 */
struct CellDescriptor
{
  /// region index
  localIndex region;
  /// subregion index
  localIndex subRegion;
  /// cell index
  localIndex index;

  /**
   * @brief Constructor for the CellDescriptor struct
   * @param[r] region index
   * @param[sr] subregion index
   * @param[i] cell index
   */
  CellDescriptor( localIndex r, localIndex sr, localIndex i )
    : region( r ), subRegion( sr ), index( i )
  {}

  /**
   * @brief Comparison operator between two CellDescriptors.
   * @param[in] other the CellDescriptor to compare with
   * @return true if they represent the same mesh element
   */
  bool operator==( CellDescriptor const & other ) const
  {
    return( region==other.region && subRegion==other.subRegion && index==other.index );
  }
};


class ProjectionEDFMHelper
{
public:

  /*
   * @brief Constructor
   */
  ProjectionEDFMHelper( MeshLevel const & mesh,
                        CellElementStencilTPFA & stencil,
                        EmbeddedSurfaceToCellStencil & edfmStencil,
                        string const & embeddedSurfaceRegionName );

  /*
   * @brief add Fracture-matrix connections to the edfmStencil
   * and remove the appropriate connections from the cellStencil
   */
  void addNonNeighboringConnections() const;

private:

  /*
   * @brief select cell faces that will host non-neighboring fracture-matrix connections
   * @param[in] subRegionFaces faces of the subRegion
   * @param[in] hostCellID id of the hostCell (tre fractured one)
   * @param[in] fracElement the index of the fracture element
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   * @return a list of the faces that need to be disconnected.
   */
  std::vector< localIndex > selectFaces( FixedOneToManyRelation const & subRegionFaces,
                                         CellDescriptor const & hostCellID,
                                         localIndex const fracElement,
                                         EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  /*
   * @brief check the intersection a fracture element and an edge
   * @param[in] fracCenter the coordinates of the center of the fracture cell
   * @param[in] fracNomral the normal vector of the fracture element
   * @param[in] edgeIdx the edge index
   * @return whether the fracture element intersects the edge or not.
   */
  bool intersection( real64 const ( &fracCenter )[3],
                     arraySlice1d< real64 const > const & fracNormal,
                     localIndex const edgeIdx ) const;

  /*
   *
   * @brief returns true is the face has only one neighbor
   * @param[in] faceIx index of the face
   * @return whether the face is a boundary face or not.
   */
  bool isBoundaryFace( localIndex const faceIdx ) const;

  // check if the center of a face is on the same side of the fracture as cell center
  bool onLargerSide( localIndex const faceIdx,
                     real64 const signedDistanceCellCenterToFrac,
                     real64 const ( &fracCenter )[3],
                     arraySlice1d< real64 const > const & fracNormal ) const;

  /*
   * @brief compute the signed distance between the fracture and as cell center
   * @param[in] hostCellID id of the host cell
   * @param[in] fracNormal the normal vector of the fracture element
   * @param[in] fracCenter the coordinates of the center of the fracture cell
   * @param[out] cellCenterToFracCenter distance vector between the cell and the fracture.
   * return The signed distance between the fracture and the cell (based on the frac normal direction)
   */
  real64 getSignedDistanceCellCenterToFracPlane( CellDescriptor const & hostCellID,
                                                 arraySlice1d< real64 const > const & fracNormal,
                                                 real64 const (&fracCenter)[3],
                                                 real64 ( &cellCenterToFracCenter )[3] ) const;

  /*
   * @brief returns true if the signed distance from the neighbor center to the frac is of the same sign as signedDistanceCellCenterToFrac
   *(computed in the host cell)
   * @@param[in] faceIx the face index
   * @param[in] signedDistanceCellCenterToFrac signed distance between the fracture and the cell (based on the frac normal direction)
   * @param[in] hostCellID id of the host cell (region, subregion and element index)
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   * return whether the neighboring cell has the fracture on the same side of the face or not.
   */
  bool neighborOnSameSide( localIndex const faceIdx,
                           real64 const signedDistanceCellCenterToFrac,
                           CellDescriptor const & hostCellID,
                           EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  /*
   * @brief given a face and its neighboring cell, return the id of the other cell
   * @param[in] faceIdx face index
   * @param[in] hostCellID id of the host cell
   * return
   */
  CellDescriptor otherCell( localIndex const faceIdx, CellDescriptor const & hostCellID ) const;

  /*
   * @brief compute the absolute transmissibility for non-neighboring F-M connection
   * @param[in] neighborCell the neighboring cell
   * @param[in] fracElement the index of the embeddedSurface elmement (the fracture element)
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   * @param[in] faceIdx the face index
   * @param[out] trans The geometric transmissibility between fracture and matrix cells.
   */
  void fractureMatrixTransmissilibility( CellDescriptor const & neighborCell,
                                         localIndex const fracElement,
                                         EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                         localIndex const faceIdx,
                                         real64 ( &trans )[2] ) const;

  /*
   * @brief add non-neighboring F-M connection with given transmissibility tothe cell stencil
   * @param[in] fracElement fracture element index
   * @param[in] cell id of the cell
   * @param[in] transmissibility geometric transmissiblity
   * @param[in] fractureSubRegion the embeddedSurfaceSubRegion
   */
  void addNonNeighboringConnection( localIndex const fracElement,
                                    CellDescriptor const & cell,
                                    real64 const (&transmissibility)[2],
                                    EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  // Private variables
  ElementRegionManager const & m_elementManager;
  ///
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodesCoord;
  ///
  arrayView2d< localIndex const > const m_edgeToNodes;
  ///
  arrayView2d< localIndex const > const m_faceToRegions;
  ///
  arrayView2d< localIndex const > const m_faceToSubRegions;
  ///
  arrayView2d< localIndex const > const m_faceToCells;
  ///
  ArrayOfArraysView< localIndex const > m_faceToEdges;
  ///
  ArrayOfArraysView< localIndex const > const m_faceToNodes;
  ///
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const m_cellCenters;
  ///
  CellElementStencilTPFA & m_cellStencil;
  ///
  EmbeddedSurfaceToCellStencil & m_edfmStencil;
  ///
  string const m_embeddedSurfaceRegionName;

};

}  // end namespace geosx


#endif // __PROJECTIONEDFMHELPER_H_
