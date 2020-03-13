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
 *
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

/**
 * @class PerforationData
 *
 * This class keeps track of all the local perforations on this rank
 */
class PerforationData : public ObjectManagerBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  explicit PerforationData( string const & name,
                            dataRepository::Group * const parent );

  /**
   * @brief default destructor
   */
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

  static string CatalogName() { return "PerforationData"; }

  virtual const string getCatalogName() const override { return CatalogName(); }

  /**
   * @brief Setter for the global number of perforations (used for well initialization)
   * @param nPerfs the global number of perforations (obtained for InternalWellGenerator)
   */
  void SetNumPerforationsGlobal( globalIndex nPerfs ) { m_numPerforationsGlobal = nPerfs; }

  /**
   * @brief Getter for the global number of perforations (used for well initialization)
   * @return the global number of perforations
   */
  globalIndex GetNumPerforationsGlobal() const { return m_numPerforationsGlobal; }

  /**
   * @brief Getter for perforation to mesh element connectivity
   * @return list of element region/subregion/index conencted to each perforation
   */
  ToElementRelation< array1d< localIndex > > & GetMeshElements() { return m_toMeshElements; }

  /**
   * @brief Const getter for perforation to mesh element connectivity
   * @return list of element region/subregion/index connected to each perforation
   */
  ToElementRelation< array1d< localIndex > > const & GetMeshElements() const { return m_toMeshElements; }

  /**
   * @brief Getter for perforation to well element connectivity
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex > & GetWellElements() { return m_wellElementIndex; }

  /**
   * @brief Const getter for perforation to well element connectivity
   * @return list of well element index connected to each perforation
   */
  arrayView1d< localIndex const > const & GetWellElements() const { return m_wellElementIndex; }

  /**
   * @brief Getter for perforation locations
   * @return list of perforation locations
   */
  arrayView1d< R1Tensor > & GetLocation() { return m_location; }

  /**
   * @brief Const getter for perforation locations
   * @return list of perforation locations
   */
  arrayView1d< R1Tensor const > const & GetLocation() const { return m_location; }

  /**
   * @brief Getter for perforation transmissibilities
   * @return list of perforation transmissibilities
   */
  arrayView1d< real64 > & GetTransmissibility() { return m_transmissibility; }

  /**
   * @brief Getter for perforation transmissibilities
   * @return list of perforation transmissibilities
   */
  arrayView1d< real64 const > const & getTransmissibility() const { return m_transmissibility; }


  /**
   * @brief Locates connected local mesh elements and resizes current object appropriately
   * @param[in] mesh the target mesh level
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   */
  void ConnectToMeshElements( MeshLevel const & mesh,
                              InternalWellGenerator const & wellGeometry );

  /**
   * @brief Connect each perforation to a local wellbore element
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] wellElementGlobalToLocalMap the global to local map of wellbore elements
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
    static constexpr auto transmissibilityString          = "transmissibility";

    dataRepository::ViewKey numPerforationsGlobal     = { numPerforationsGlobalString };
    dataRepository::ViewKey reservoirElementRegion    = { reservoirElementRegionString };
    dataRepository::ViewKey reservoirElementSubregion = { reservoirElementSubregionString };
    dataRepository::ViewKey reservoirElementIndex     = { reservoirElementIndexString };
    dataRepository::ViewKey wellElementIndex          = { wellElementIndexString };
    dataRepository::ViewKey location                  = { locationString };
    dataRepository::ViewKey transmissibility          = { transmissibilityString };

  } viewKeysPerforationData;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {} groupKeysPerforationData;

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const group ) override;

private:

  void DebugLocalPerforations() const;

  /// global number of perforations
  globalIndex m_numPerforationsGlobal;

  /// indices of the mesh elements connected to perforations
  ToElementRelation< array1d< localIndex > > m_toMeshElements;

  /// indices of the well element to which perforations are attached
  array1d< localIndex > m_wellElementIndex;

  /// location of the perforations
  array1d< R1Tensor > m_location;

  /// transmissibility (well index) of the perforations
  array1d< real64 > m_transmissibility;

};

} //namespace geosx

#endif //GEOSX_WELLS_PERFORATIONDATA_HPP
