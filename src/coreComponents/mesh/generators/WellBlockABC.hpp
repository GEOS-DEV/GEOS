/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_WELLBLOCKABC_HPP
#define GEOSX_WELLBLOCKABC_HPP

#include "dataRepository/Group.hpp"
#include "mesh/ElementType.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"

#include <vector>

namespace geosx
{

/**
 * Abstract base class defining the information provided by any well block implementation.
 *
 * It's noteworthy that the WellblockABC is immutable oriented.
 * The derived implementations need to have the modification/creation capabilities.
 */
class WellBlockABC
{
public:

  /**
   * @name Getters / Setters
   */
  ///@{

  // getters for element data

  /**
   * @brief Get the global number of well elements.
   * @return the global number of elements
   */
  virtual globalIndex numElements() const = 0;

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  virtual arrayView2d< real64 const > elemCoords() const = 0;

  /**
   * @brief Get the global indices mapping an element to the next.
   * @return list providing the global index of the next element for each element
   */
  virtual arrayView1d< globalIndex const > nextElemIndex() const = 0;

  /**
   * @brief Get the global indices mapping an element to the previous ones.
   * @return list providing the global indices of the previous elements for each element
   */
  virtual arrayView1d< arrayView1d< globalIndex const > const > prevElemIndices() const = 0;

  /**
   * @brief Get the global indices of the well nodes nodes connected to each element.
   * @return list providing the global index of the well nodes for each well element
   */
  virtual arrayView2d< globalIndex const > elemToNodesMap() const = 0;

  /**
   * @brief Get the volume of the well elements.
   * @return list of volumes of the well elements
   */
  virtual arrayView1d< real64 const > elemVolume() const = 0;

  /**
   * @brief Get the radius in the well.
   * @return the radius in the well
   */
  virtual real64 elementRadius() const = 0;

  // getters for node data

  /**
   * @brief Get the global number of well nodes.
   * @return the global number of nodes
   */
  virtual globalIndex numNodes() const = 0;

  /**
   * @brief Get the physical location of the centers of well elements.
   * @return list of center locations of the well elements
   */
  virtual arrayView2d< real64 const > nodeCoords() const = 0;



  // getters for perforation data

  /**
   * @brief Get the global number of perforations on this well.
   * @return the global number of elements
   */
  virtual globalIndex numPerforations() const = 0;

  /**
   * @brief Get the locations of the perforations.
   * @return list of locations of all the perforations on the well
   */
  virtual arrayView2d< real64 const > perfCoords() const = 0;

  /**
   * @brief Get the well transmissibility at the perforations.
   * @return list of well transmissibility at all the perforations on the well
   */
  virtual arrayView1d< real64 const > perfTransmissibility() const = 0;

  /**
   * @brief Get the global indices of the well elements connected to each perforation.
   * @return list providing the global index of the connected well element for each perforation
   */
  virtual arrayView1d< globalIndex const > perfElemIndex() const = 0;

};

}

#endif //GEOSX_WELLBLOCKABC_HPP
