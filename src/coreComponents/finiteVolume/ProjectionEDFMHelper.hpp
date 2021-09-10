#ifndef __PROJECTIONEDFMHELPER_H_
#define __PROJECTIONEDFMHELPER_H_

#include <list>                                  // provides std::vector
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
  bool operator==( CellDescriptor const & other )
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

  virtual ~ProjectionEDFMHelper() = default;

private:

  /*
   * @brief select cell faces that will host non-neighboring fracture-matrix connections
   */
  std::list< localIndex > selectFaces( FixedOneToManyRelation const & subRegionFaces,
                                       CellDescriptor const & hostCellID,
                                       localIndex const fracElement,
                                       EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  /*
   * @brief check the intersection  a fracture element and an edge
   */
  bool intersection( real64 const ( &fracOrigin )[3],
                     arraySlice1d< real64 const > const & fracNormal,
                     localIndex const edgeIdx,
                     real64 ( &tmp )[3] ) const;

  /*
   *
   * @brief returns true is the face has only one neighbor
   */
  bool isBoundaryFace( localIndex const faceIdx ) const;

  // check if the center of a face is on the same side of the fracture as cell center
  bool onLargerSide( localIndex const faceIdx,
                     real64 const signedDistanceCellCenterToFrac,
                     real64 const ( &fracOrigin )[3],
                     arraySlice1d< real64 const > const & fracNormal ) const;

  /*
   * @brief compute the signed distance between the fracture and as cell center
   */
  real64 getSignedDistanceCellCenterToFracPlane( CellDescriptor const & hostCellID,
                                                 arraySlice1d< real64 const > const & fracNormal,
                                                 real64 const (&fracOrigin)[3],
                                                 real64 ( &tmp )[3] ) const;

  /*
   * @brief returns true if the signed distance from the neighbor center to the frac is of the same sign as signedDistanceCellCenterToFrac
   *(computed in the host cell)
   */
  bool neighborOnSameSide( localIndex const faceIdx,
                           real64 const signedDistanceCellCenterToFrac,
                           CellDescriptor const & hostCellID,
                           EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  /*
   * @brief given a face and its neighboring cell, return the id of the other cell
   */
  CellDescriptor otherCell( localIndex const faceIdx, CellDescriptor const & hostCellID ) const;

  /*
   * @brief compute the absolute transmissibility for non-neighboring F-M connection
   */
  void fractureMatrixTransmissilibility( CellDescriptor const & neighborCell,
                                         localIndex const fracElement,
                                         EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                         localIndex const faceIdx,
                                         real64 (& trans)[2] ) const;

  /*
   * @brief add non-neighboring F-M connection with given transmissibility tothe cell stencil
   */
  void addNonNeighboringConnection( localIndex const fracElement,
                                    CellDescriptor const & cell,
                                    real64 const (& transmissibility)[2],
                                    EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  // Private variables
  ElementRegionManager const & m_elementManager;
  ///
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodesCoord;
  ///
  arrayView2d< localIndex const > const m_edgeToNodes;
  ///
  arrayView2d< localIndex const > const m_facesToCells;
  ///
  arrayView2d< localIndex const > const m_facesToRegions;
  ///
  arrayView2d< localIndex const > const m_facesToSubRegions;
  ///
  ArrayOfArraysView< localIndex const > const m_facesToNodes;
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
