/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file CapillaryPressureBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CAPILLARYPRESSUREBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CAPILLARYPRESSUREBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

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
			 dataRepository::ManagedGroup * const parent );
  
  virtual ~CapillaryPressureBase() override;

  // *** ManagedGroup interface
  
  virtual void ProcessInputFile_PostProcess() override;
  
  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** CapillaryPressure-specific interface
  
  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] phaseVolFraction input phase volume fraction
   */
  virtual void BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction ) = 0;

  /**
   * @brief Perform a single point constitutive update.
   * @param[in] phaseVolFraction input phase volume fraction
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] q second constitutive index (e.g. quadrature index)
   *
   * @note This function should generally not be called from a kernel, use BatchUpdate instead
   */
  virtual void PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) {}

  localIndex numFluidPhases() const;

  string const & phaseName( localIndex ip ) const;

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phaseNamesString     = "phaseNames";
    static constexpr auto phaseTypesString     = "phaseTypes";
    static constexpr auto phaseOrderString     = "phaseOrder";
    
    static constexpr auto phaseCapPressureString                    = "phaseCapPressure";                    // Pc_p
    static constexpr auto dPhaseCapPressure_dPhaseVolFractionString = "dPhaseCapPressure_dPhaseVolFraction"; // dPc_p/dS_p

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseNames = { phaseNamesString };
    ViewKey phaseTypes = { phaseTypesString };
    ViewKey phaseOrder = { phaseOrderString };

    ViewKey phaseCapPressure                    = { phaseCapPressureString };                    // Pc_p
    ViewKey dPhaseCapPressure_dPhaseVolFraction = { dPhaseCapPressure_dPhaseVolFractionString }; // dPc_p/dS_p

  } viewKeysCapillaryPressureBase;

protected:

  /**
   * @brief Function to batch process constitutive updates via a kernel launch.
   * @tparam LEAFCLASS The derived class that provides the functions for use in the kernel
   * @tparam ARGS Parameter pack for arbitrary number of arbitrary types for the function parameter list
   * @param phaseVolumeFraction array containing the phase volume fraction, which is input to the update.
   * @param args arbitrary number of arbitrary types that are passed to the kernel
   */
  template< typename LEAFCLASS, typename ... ARGS >
  void BatchUpdateKernel( arrayView2d<real64 const> const & phaseVolumeFraction,
                          ARGS&& ... args );

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void ResizeFields( localIndex const size, localIndex const numPts );

  // phase names read from input
  string_array     m_phaseNames;

  // phase ordering info
  array1d<integer> m_phaseTypes;
  array1d<integer> m_phaseOrder;
  
  // output quantities
  array3d<real64>  m_phaseCapPressure;
  array4d<real64>  m_dPhaseCapPressure_dPhaseVolFrac;
  
};

template< typename LEAFCLASS, typename ... ARGS >
void CapillaryPressureBase::BatchUpdateKernel( arrayView2d<real64 const> const & phaseVolumeFraction,
                                               ARGS&& ... args )
{
  localIndex const numElem = m_phaseCapPressure.size(0);
  localIndex const numQ    = m_phaseCapPressure.size(1);
  localIndex const NP      = numFluidPhases();

  arrayView3d<real64> const & phaseCapPressure = m_phaseCapPressure;
  arrayView4d<real64> const & dPhaseCapPressure_dPhaseVolFrac = m_dPhaseCapPressure_dPhaseVolFrac;
  arrayView1d<integer const> const & phaseOrder = m_phaseOrder;
  
  forall_in_range( 0, numElem, GEOSX_LAMBDA ( localIndex const k )
  {
    for( localIndex q=0 ; q<numQ ; ++q )
    {
      LEAFCLASS::Compute( NP,
			  phaseVolumeFraction[k],
			  phaseCapPressure[k][q],
			  dPhaseCapPressure_dPhaseVolFrac[k][q],
			  phaseOrder,
			  args... );
    }
  });
}

  
} // namespace constitutive

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CAPILLARYPRESSUREBASE_HPP
