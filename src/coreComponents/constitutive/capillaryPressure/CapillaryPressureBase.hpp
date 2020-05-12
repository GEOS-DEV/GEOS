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

/**
 * @file CapillaryPressureBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREBASE_HPP
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

class CapillaryPressureBase : public ConstitutiveBase
{
public:

  struct PhaseType
  {
    static constexpr integer OIL            = 0;
    static constexpr integer GAS            = 1;
    static constexpr integer WATER          = 2;
    static constexpr integer MAX_NUM_PHASES = 3;
  };

  // choose the reference pressure to be the oil pressure for all models
  static constexpr integer REFERENCE_PHASE = PhaseType::OIL;

  CapillaryPressureBase( std::string const & name,
                         dataRepository::Group * const parent );

  virtual ~CapillaryPressureBase() override;

  // *** Group interface

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** CapillaryPressure-specific interface

  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] phaseVolFraction input phase volume fraction
   */
  virtual void BatchUpdate( arrayView2d< real64 const > const & phaseVolumeFraction ) = 0;

  /**
   * @brief Perform a single point constitutive update.
   * @param[in] phaseVolFraction input phase volume fraction
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] q second constitutive index (e.g. quadrature index)
   *
   * @note This function should generally not be called from a kernel, use BatchUpdate instead
   */
  virtual void PointUpdate( arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( phaseVolFraction ),
                            localIndex const GEOSX_UNUSED_PARAM( k ),
                            localIndex const GEOSX_UNUSED_PARAM( q ) ) {}

  localIndex numFluidPhases() const { return m_phaseNames.size(); }

  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  arrayView3d< real64 const > phaseCapPressure() const { return m_phaseCapPressure; }
  arrayView4d< real64 const > dPhaseCapPressure_dPhaseVolFraction() const { return m_dPhaseCapPressure_dPhaseVolFrac; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phaseNamesString = "phaseNames";
    static constexpr auto phaseTypesString = "phaseTypes";
    static constexpr auto phaseOrderString = "phaseOrder";

    static constexpr auto phaseCapPressureString                    = "phaseCapPressure";                    // Pc_p
    static constexpr auto dPhaseCapPressure_dPhaseVolFractionString = "dPhaseCapPressure_dPhaseVolFraction"; // dPc_p/dS_p
  } viewKeysCapillaryPressureBase;

protected:
  virtual void PostProcessInput() override;

  /**
   * @brief Function to batch process constitutive updates via a kernel launch.
   * @tparam LEAFCLASS The derived class that provides the functions for use in the kernel
   * @tparam ARGS Parameter pack for arbitrary number of arbitrary types for the function parameter list
   * @param phaseVolumeFraction array containing the phase volume fraction, which is input to the update.
   * @param args arbitrary number of arbitrary types that are passed to the kernel
   */
  template< typename LEAFCLASS, typename POLICY=serialPolicy, typename ... ARGS >
  void BatchUpdateKernel( arrayView2d< real64 const > const & phaseVolumeFraction,
                          ARGS && ... args );

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void ResizeFields( localIndex const size, localIndex const numPts );

  // phase names read from input
  string_array m_phaseNames;

  // phase ordering info
  array1d< integer > m_phaseTypes;
  array1d< integer > m_phaseOrder;

  // output quantities
  array3d< real64 >  m_phaseCapPressure;
  array4d< real64 >  m_dPhaseCapPressure_dPhaseVolFrac;

};

template< typename LEAFCLASS, typename POLICY, typename ... ARGS >
void CapillaryPressureBase::BatchUpdateKernel( arrayView2d< real64 const > const & phaseVolumeFraction,
                                               ARGS && ... args )
{
  localIndex const numElem = m_phaseCapPressure.size( 0 );
  localIndex const numQ    = m_phaseCapPressure.size( 1 );
  localIndex const NP      = numFluidPhases();

  arrayView3d< real64 > const & phaseCapPressure = m_phaseCapPressure;
  arrayView4d< real64 > const & dPhaseCapPressure_dPhaseVolFrac = m_dPhaseCapPressure_dPhaseVolFrac;
  arrayView1d< integer const > const & phaseOrder = m_phaseOrder;

  forAll< POLICY >( numElem, [=] ( localIndex const k )
  {
    for( localIndex q=0; q<numQ; ++q )
    {
      LEAFCLASS::Compute( NP,
                          phaseVolumeFraction[k],
                          phaseCapPressure[k][q],
                          dPhaseCapPressure_dPhaseVolFrac[k][q],
                          phaseOrder,
                          args ... );
    }
  } );
}


} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREBASE_HPP
