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
 * @file RelativePermeabilityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

class RelativePermeabilityBase : public ConstitutiveBase
{
public:

  struct PhaseType
  {
    static constexpr integer OIL            = 0;
    static constexpr integer GAS            = 1;
    static constexpr integer WATER          = 2;
    static constexpr integer MAX_NUM_PHASES = 3;
  };

  // order of the phase properties in the water-oil data
  struct WaterOilPairPhaseType
  {
    static constexpr integer WATER = 0; // first water phase property
    static constexpr integer OIL   = 1; // second oil phase property
  };

  // order of the phase properties in the gas-oil data
  struct GasOilPairPhaseType
  {
    static constexpr integer GAS   = 0; // first gas phase property
    static constexpr integer OIL   = 1; // second oil phase property
  };



  RelativePermeabilityBase( std::string const & name, dataRepository::Group * const parent );

  virtual ~RelativePermeabilityBase() override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Function to update state of a single material point.
   * @param[in] phaseVolFraction input phase volume fraction
   * @param[in] k the first index of the storage arrays (elem index)
   * @param[in] q the secound index of the storage arrays (quadrature index)
   *
   * @note This function performs a point update, but should not be called
   *       within a kernel since it is virtual, and the required data is not
   *       guaranteed to be in the target memory space.
   */
  virtual void PointUpdate( arraySlice1d< real64 const > const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) = 0;

  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] phaseVolFraction input phase volume fraction
   */
  virtual void BatchUpdate( arrayView2d< real64 const > const & phaseVolumeFraction ) = 0;

  localIndex numFluidPhases() const { return m_phaseNames.size(); }

  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  arrayView3d< real64 const > phaseRelPerm() const { return m_phaseRelPerm; }
  arrayView4d< real64 const > dPhaseRelPerm_dPhaseVolFraction() const { return m_dPhaseRelPerm_dPhaseVolFrac; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phaseNamesString     = "phaseNames";
    static constexpr auto phaseTypesString     = "phaseTypes";
    static constexpr auto phaseOrderString     = "phaseTypes";

    static constexpr auto phaseRelPermString                    = "phaseRelPerm";                    // Kr
    static constexpr auto dPhaseRelPerm_dPhaseVolFractionString = "dPhaseRelPerm_dPhaseVolFraction"; // dKr_p/dS_p
  } viewKeysRelativePermeabilityBase;

protected:
  virtual void PostProcessInput() override;

  /**
   * @brief Function to batch process constitutive updates via a kernel launch.
   * @tparam LEAFCLASS The derived class that provides the functions for use
   *                   in the kernel
   * @tparam ARGS Parameter pack for arbitrary number of arbitrary types for
   *              the function parameter list
   * @param phaseVolumeFraction array containing the phase volume fraction,
   *                            which is input to the update.
   * @param args arbitrary number of arbitrary types that are passed to the
   *             kernel
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
  array3d< real64 >  m_phaseRelPerm;
  array4d< real64 >  m_dPhaseRelPerm_dPhaseVolFrac;


};


template< typename LEAFCLASS, typename POLICY, typename ... ARGS >
void RelativePermeabilityBase::BatchUpdateKernel( arrayView2d< real64 const > const & phaseVolumeFraction,
                                                  ARGS && ... args )
{
  localIndex const numElem = m_phaseRelPerm.size( 0 );
  localIndex const numQ = m_phaseRelPerm.size( 1 );
  localIndex const NP = numFluidPhases();

  arrayView3d< real64 > const & phaseRelPerm = m_phaseRelPerm;
  arrayView4d< real64 > const & dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac;

  forAll< POLICY >( numElem, [=] ( localIndex const k )
  {
    for( localIndex q=0; q<numQ; ++q )
    {
      LEAFCLASS::Compute( NP,
                          phaseVolumeFraction[k],
                          phaseRelPerm[k][q],
                          dPhaseRelPerm_dPhaseVolFrac[k][q],
                          args ... );
    }
  } );
}

} // namespace constitutive

} // namespace geosx


#endif //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
