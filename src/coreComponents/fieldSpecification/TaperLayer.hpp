/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * TaperLayer.hpp
 *
 */

#ifndef GEOS_FIELDSPECIFICATION_TAPERLAYER_HPP_
#define GEOS_FIELDSPECIFICATION_TAPERLAYER_HPP_

#include "FieldSpecificationBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

class DomainPartition;

/**
 * @class TaperLayer
 * A class to manage Perfectly Matched Layer for wave propagation solvers
 */
class TaperLayer : public FieldSpecificationBase
{
public:
  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationBase in the data repository
   * @param parent the parent group of this group.
   */
  TaperLayer( string const & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  TaperLayer() = delete;

  /**
   * @brief destructor
   */
  virtual ~TaperLayer() = default;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "Taper"; }

  /**
   * @brief Getter for the taper minimum coordinates
   * @return Mininum (x,y,z) coordinates of inner taper boundaries
   */
  R1Tensor32 getMin() const { return m_xMin; }

  /**
   * @brief Getter for the taper maximum coordinates
   * @return Maximum (x,y,z) coordinates of inner taper boundaries
   */
  R1Tensor32 getMax() const { return m_xMax; }

  /**
   * @brief Getter for the taper thickness
   * @return Thickness of the taper region at left, front, and top sides
   */
  R1Tensor32 getThicknessMinXYZ() const { return m_thicknessMinXYZ; }

  /**
   * @brief Getter for the taper thickness
   * @return Thickness of the taper region at right, back, and bottom sides
   */
  R1Tensor32 getThicknessMaxXYZ() const { return m_thicknessMaxXYZ; }

  /**
   *
   */
  real32 getTaperCoeff() const { return m_taperCoeff; }

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationBase::viewKeyStruct
  {
    /// @return String key for the mininum (x,y,z) coordinates of inner taper boundaries
    static constexpr char const * xMinString() { return "xMin"; }

    /// @return String key for the maximum (x,y,z) coordinates of inner taper boundaries
    static constexpr char const * xMaxString() { return "xMax"; }

    /// @return String key for the reflectivity of the PML region
    static constexpr char const * taperCoeffString() { return "taperCoeff"; }

    /// @return String key for the thickness of the PML region (left, front, top sides)
    static constexpr char const * thicknessMinXYZString() { return "thicknessMinXYZ"; }

    /// @return String key for the thickness of the PML region (right, back, bottom sides)
    static constexpr char const * thicknessMaxXYZString() { return "thicknessMaxXYZ"; }


  };

  /**
   * @brief Safeguard for the minimum allowed taper thickness
   */
  static constexpr real64 minThickness = 0.001;

  /**
   * @brief Smallest possible values for xMin, below which they are computed internally
   */
  static constexpr real64 smallestXMin = -999999.0;

  /**
   * @brief Largest possible values for xMax, below which they are computed internally
   */
  static constexpr real64 largestXMax = 999999.0;

protected:

  virtual void postProcessInput() override final;

private:

  /// Mininum (x,y,z) coordinates of inner taper boundaries
  R1Tensor32 m_xMin;

  /// Maximum (x,y,z) coordinates of inner taper boundaries
  R1Tensor32 m_xMax;

  /// Thickness of the taper region, used to compute the damping profile
  R1Tensor32 m_thicknessMinXYZ;
  R1Tensor32 m_thicknessMaxXYZ;

  // Coefficient for the taper, used to compile damping profile
  real32 m_taperCoeff;

};

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_TAPER_LAYER_HPP_ */
