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
 * PerfectlyMatchedLayer.hpp
 *
 */

#ifndef GEOSX_FIELDSPECIFICATION_PERFECTLYMATCHEDLAYER_HPP_
#define GEOSX_FIELDSPECIFICATION_PERFECTLYMATCHEDLAYER_HPP_

#include "FieldSpecificationBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geosx
{

class DomainPartition;

/**
 * @class PerfectlyMatchedLayer
 * A class to manage Perfectly Matched Layer for wave propagation solvers
 */
class PerfectlyMatchedLayer : public FieldSpecificationBase
{
public:
  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationBase in the data repository
   * @param parent the parent group of this group.
   */
  PerfectlyMatchedLayer( string const & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  PerfectlyMatchedLayer() = delete;

  /**
   * @brief destructor
   */
  virtual ~PerfectlyMatchedLayer() = default;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "PML"; }

  /**
   * @brief Getter for the PML minimum coordinates
   * @return Mininum (x,y,z) coordinates of inner PML boundaries
   */
  R1Tensor getMin() const { return m_xMin; }

  /**
   * @brief Copy for the PML minimum coordinates
   * @param arr Mininum (x,y,z) coordinates of inner PML boundaries
   */
  void copyMin( real64 (& arr)[3] ) const { for( int i=0; i<3; ++i ) arr[i]=m_xMin[i]; }

  /**
   * @brief Getter for the PML maximum coordinates
   * @return Maximum (x,y,z) coordinates of inner PML boundaries
   */
  R1Tensor getMax() const { return m_xMax; }

  /**
   * @brief Copy for the PML maximum coordinates
   * @param arr Maximum (x,y,z) coordinates of inner PML boundaries
   */
  void copyMax( real64 (& arr)[3] ) const { for( int i=0; i<3; ++i ) arr[i]=m_xMax[i]; }

  /**
   * @brief Getter for the PML reflectivity
   * @return Desired reflectivity of the PML region
   */
  real64 getReflectivity() const { return m_reflectivity; }

  /**
   * @brief Getter for the PML thickness
   * @return Thickness of the PML region
   */
  R1Tensor getThicknessMinXYZ() const { return m_thicknessMinXYZ; }
  R1Tensor getThicknessMaxXYZ() const { return m_thicknessMaxXYZ; }

  /**
   * @brief Copy for the PML thickness
   * @param arr Thickness of the PML region
   */
  void copyThicknessMinXYZ( real64 (& arr)[3] ) const { for( int i=0; i<3; ++i ) arr[i]=m_thicknessMinXYZ[i]; }
  void copyThicknessMaxXYZ( real64 (& arr)[3] ) const { for( int i=0; i<3; ++i ) arr[i]=m_thicknessMaxXYZ[i]; }

  /**
   * @brief Getter for the PML wave speed
   * @return Wave speed of the PML region
   */
  R1Tensor getWaveSpeedMinXYZ() const { return m_waveSpeedMinXYZ; }
  R1Tensor getWaveSpeedMaxXYZ() const { return m_waveSpeedMaxXYZ; }

  /**
   * @brief Copy for the PML wave speed
   * @param arr Wave speed of the PML region
   */
  void copyWaveSpeedMinXYZ( real64 (& arr)[3] ) const { for( int i=0; i<3; ++i ) arr[i]=m_waveSpeedMinXYZ[i]; }
  void copyWaveSpeedMaxXYZ( real64 (& arr)[3] ) const { for( int i=0; i<3; ++i ) arr[i]=m_waveSpeedMaxXYZ[i]; }

  /**
   * @brief Setter for the PML minimum coordinates
   * @param arr Mininum (x,y,z) coordinates of inner PML boundaries
   */
  void setMin( real64 const (&arr)[3] ) { for( int i=0; i<3; ++i ) m_xMin[i]=arr[i]; }

  /**
   * @brief Setter for the PML maximum coordinates
   * @param arr Maximum (x,y,z) coordinates of inner PML boundaries
   */
  void setMax( real64 const (&arr)[3] ) { for( int i=0; i<3; ++i ) m_xMax[i]=arr[i]; }

  /**
   * @brief Setter for the PML thickness
   * @param arr Thickness of the PML region
   */
  void setThicknessMinXYZ( real64 const (&arr)[3] ) { for( int i=0; i<3; ++i ) m_thicknessMinXYZ[i]=arr[i]; }
  void setThicknessMaxXYZ( real64 const (&arr)[3] ) { for( int i=0; i<3; ++i ) m_thicknessMaxXYZ[i]=arr[i]; }

  /**
   * @brief Setter for the PML wave speed
   * @param arr Wave speed of the PML region
   */
  void setWaveSpeedMinXYZ( real64 arr[3] ) { for( int i=0; i<3; ++i ) m_waveSpeedMinXYZ[i]=arr[i]; }
  void setWaveSpeedMaxXYZ( real64 arr[3] ) { for( int i=0; i<3; ++i ) m_waveSpeedMaxXYZ[i]=arr[i]; }

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationBase::viewKeyStruct
  {
    static constexpr char const * xMinString() { return "xMin"; }
    static constexpr char const * xMaxString() { return "xMax"; }
    static constexpr char const * reflectivityString() { return "reflectivity"; }
    static constexpr char const * thicknessMinXYZString() { return "thicknessMinXYZ"; }
    static constexpr char const * thicknessMaxXYZString() { return "thicknessMaxXYZ"; }
    static constexpr char const * waveSpeedMinXYZString() { return "waveSpeedMinXYZ"; }
    static constexpr char const * waveSpeedMaxXYZString() { return "waveSpeedMaxXYZ"; }

  };

protected:

  virtual void postProcessInput() override final;

  virtual void initializePreSubGroups() override final;

private:

  /// Mininum (x,y,z) coordinates of inner PML boundaries
  R1Tensor m_xMin;

  /// Maximum (x,y,z) coordinates of inner PML boundaries
  R1Tensor m_xMax;

  /// Desired reflectivity of the PML region, used to compute the damping profile
  real64 m_reflectivity;

  /// Thickness of the PML region, used to compute the damping profile
  R1Tensor m_thicknessMinXYZ;
  R1Tensor m_thicknessMaxXYZ;

  /// Wave speed in the PML region, used to compute the damping profile
  R1Tensor m_waveSpeedMinXYZ;
  R1Tensor m_waveSpeedMaxXYZ;

};

} /* namespace geosx */

#endif /* GEOSX_FIELDSPECIFICATION_PERFECTLYMATCHEDLAYER_HPP_ */
