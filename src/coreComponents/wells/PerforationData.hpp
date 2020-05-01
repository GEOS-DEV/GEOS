/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file PerforationData.hpp
 */

#ifndef GEOSX_WELLS_PERFORATIONDATA_HPP
#define GEOSX_WELLS_PERFORATIONDATA_HPP

#include "dataRepository/Group.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "mesh/ToElementRelation.hpp"
#include "InternalWellGenerator.hpp"

namespace geosx
{

class DomainPartition;
class MeshLevel;
class WellElementSubRegion;
class CellBlock;

/**
 * @class PerforationData
 *
 * This class keeps track of all the local perforations on this rank
 */
class PerforationData : public ObjectManagerBase
{
public:

  /**
   * @brief main constructor for PerforationData Objects
   * @param [in] name the name of this instantiation of PerforationData in the repository
   * @param [in] parent the parent group of this instantiation of PerforationData
   */
  explicit PerforationData( string const & name, dataRepository::Group * const parent );


  /// default destructor
  ~PerforationData() override;


  /// deleted default constructor
  PerforationData() = delete;


  /// deleted copy constructor
  PerforationData( PerforationData const & ) = delete;

  
  /// deleted move constructor
  PerforationData( PerforationData && ) = delete;

  
  /// deleted assignment operator
  PerforationData & operator=( PerforationData const & ) = delete;

  
  /// deleted move operator
  PerforationData & operator=( PerforationData && ) = delete;

  
  /// Get the catalog name
  static string CatalogName() { return "PerforationData"; }

  
  /// Get a const catalog name
  virtual const string getCatalogName() const override { return CatalogName(); }

  
  /**
   * @brief Set the global number of perforations used for well initialization.
   * @param [in] nPerfs global number of perforations (obtained for InternalWellGenerator)
   */
  void SetNumPerforationsGlobal( globalIndex nPerfs ) { m_numPerforationsGlobal = nPerfs; }

  
  /**
   * @brief Get the global number of perforations (used for well initialization).
   * @return global number of perforations
   */
  globalIndex GetNumPerforationsGlobal() const { return m_numPerforationsGlobal; }

  
  /**
   * @brief Get perforation-to-mesh-element connectivity.
   * @return list of element region/subregion/index connected to each perforation
   */
  ToElementRelation< array1d< localIndex > > & GetMeshElements() { return m_toMeshElements; }

  
  /**
   * @brief Provides an immutable accessor to a const perforation-to-mesh-element connectivity.
   * @return list of element region/subregion/index connected to each perforation
   */
  ToElementRelation< array1d< localIndex > > const & GetMeshElements() const { return m_toMeshElements; }

  
  /**
   * @brief Get perforation-to-well-element connectivity.
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex > & GetWellElements() { return m_wellElementIndex; }

  
  /**
   * @brief Provides an immutable accessor to a const perforation-to-well-element connectivity.
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex const > const & GetWellElements() const { return m_wellElementIndex; }

  
  /**
   * @brief Get perforation locations.
   * @return list of perforation locations
   */
  arrayView1d< R1Tensor > & GetLocation() { return m_location; }

  
  /**
   * @brief Provides an immutable accessor to a const perforation location arrayView.
   * @return list of perforation locations
   */
  arrayView1d< R1Tensor const > const & GetLocation() const { return m_location; }

  
  /**
   * @brief Provides an immutable accessor to a const perforation well indices array.
   * @return list of perforation well indices
   */
  arrayView1d< real64 const > const & GetWellTransmissibility() const { return m_wellTransmissibility; }

  
  /**
   * @brief Get perforation well indices.
   * @return list of perforation well indices
   */
  arrayView1d< real64 > & GetWellTransmissibility() { return m_wellTransmissibility; }

  
  /**
   * @brief Computes the well transmissibility for each local perforation on this well.
   * @param[in] mesh target mesh level
   * @param[in] wellElementSubRegion  subRegion corresponding to this well
   * @param[in] permeabilityKey key to access the permeability in the reservoir
   */
  void ComputeWellTransmissibility( MeshLevel const & mesh,
                                    WellElementSubRegion const * const wellElemSubRegion,
                                    string const & permeabilityKey );

  
  /**
   * @brief Locates connected local mesh elements and resizes current object appropriately.
   * @param[in] mesh target mesh level
   * @param[in] wellGeometry  InternalWellGenerator containing the global well topology
   */
  void ConnectToMeshElements( MeshLevel const & mesh,
                              InternalWellGenerator const & wellGeometry );

  
  /**
   * @brief Connects each perforation to a local wellbore element.
   * @param[in] wellGeometry InternalWellGenerator containing the global well topology
   * @param[in] wellElementGlobalToLocalMap  global-to-local map of wellbore elements
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */
  void ConnectToWellElements( InternalWellGenerator const & wellGeometry,
                              unordered_map< globalIndex, localIndex > const & globalToLocalWellElementMap,
                              globalIndex elemOffsetGlobal );

  
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto numPerforationsGlobalString     = "numPerforationsGlobal";
    static constexpr auto reservoirElementRegionString    = "reservoirElementRegion";
    static constexpr auto reservoirElementSubregionString = "reservoirElementSubregion";
    static constexpr auto reservoirElementIndexString     = "reservoirElementIndex";
    static constexpr auto wellElementIndexString          = "wellElementIndex";
    static constexpr auto locationString                  = "location";
    static constexpr auto wellTransmissibilityString      = "wellTransmissibility";

    dataRepository::ViewKey numPerforationsGlobal     = { numPerforationsGlobalString };
    dataRepository::ViewKey reservoirElementRegion    = { reservoirElementRegionString };
    dataRepository::ViewKey reservoirElementSubregion = { reservoirElementSubregionString };
    dataRepository::ViewKey reservoirElementIndex     = { reservoirElementIndexString };
    dataRepository::ViewKey wellElementIndex          = { wellElementIndexString };
    dataRepository::ViewKey location                  = { locationString };
    dataRepository::ViewKey wellTransmissibility      = { wellTransmissibilityString };

  } viewKeysPerforationData;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {} groupKeysPerforationData;


private:

  /**
   * @brief Computes the approximate dimensions of the reservoir element containing a perforation.
   *        This is done by computing a bounding box containing the element
   * @param[in] mesh target mesh level
   * @param[in] er  index of the element region containing the reservoir element
   * @param[in] esr  index of the element subRegion containing the reservoir element
   * @param[in] ei  index of the reservoir element
   * @param[inout] dx dimension of the element in the x-direction
   * @param[inout] dy dimension of the element in the y-direction
   * @param[inout] dz dimension of the element in the z-direction
   */
  void GetReservoirElementDimensions( MeshLevel const & mesh,
                                      localIndex const er, localIndex const esr, localIndex const ei,
                                      real64 & dx, real64 & dy, real64 & dz ) const;

  /**
   * @brief Checks if the well is along the x-, y-, or z- directions.
   * @param[in] vecWellElemCenterToPerf vector connecting the well element center to the perforation
   * @param[in] dx dimension of the element in the x-direction
   * @param[in] dy dimension of the element in the y-direction
   * @param[in] dz dimension of the element in the z-direction
   * @param[in] perm absolute permeability in the reservoir element
   * @param[inout] d1 dimension of the element in the first direction
   * @param[inout] d2 dimension of the element in the second direction
   * @param[inout] h dimension of the element in the third direction
   * @param[inout] k1 absolute permeability in the reservoir element (first direction)
   * @param[inout] k2 absolute permeability in the reservoir element (second direction)
   */
  void DecideWellDirection( R1Tensor const & vecWellElemCenterToPerf,
                            real64 const & dx, real64 const & dy, real64 const & dz,
                            R1Tensor const & perm,
                            real64 & d1, real64 & d2, real64 & h,
                            real64 & k1, real64 & k2 ) const;

  /// Console output of perforations for debugging purposes
  void DebugLocalPerforations() const;

  /// global number of perforations
  globalIndex m_numPerforationsGlobal;

  /// indices of the mesh elements connected to perforations
  ToElementRelation< array1d< localIndex > > m_toMeshElements;

  /// indices of the well element to which perforations are attached
  array1d< localIndex > m_wellElementIndex;

  /// location of the perforations
  array1d< R1Tensor > m_location;

  /// well index of the perforations
  array1d< real64 > m_wellTransmissibility;

};

} //namespace geosx

#endif //GEOSX_WELLS_PERFORATIONDATA_HPP
