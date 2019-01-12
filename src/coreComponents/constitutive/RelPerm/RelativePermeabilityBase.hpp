/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
  * @file RelativePermeabilityBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

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

  RelativePermeabilityBase( std::string const & name, dataRepository::ManagedGroup * const parent );

  virtual ~RelativePermeabilityBase() override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction ) = 0;

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
  virtual void StateUpdatePointRelPerm( arraySlice1d<real64 const> const & phaseVolFraction,
                                        localIndex const k,
                                        localIndex const q ) {}

  localIndex numFluidPhases() const;

  string const & phaseName( localIndex ip ) const;

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phaseNamesString     = "phaseNames";
    static constexpr auto phaseTypesString     = "phaseTypes";
    static constexpr auto phaseOrderString     = "phaseTypes";

    static constexpr auto phaseRelPermString                    = "phaseRelPerm";                    // Kr
    static constexpr auto dPhaseRelPerm_dPhaseVolFractionString = "dPhaseRelPerm_dPhaseVolFraction"; // dKr_p/dS_p

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseNames = { phaseNamesString };
    ViewKey phaseTypes = { phaseTypesString };
    ViewKey phaseOrder = { phaseOrderString };

    ViewKey phaseRelPerm                    = { phaseRelPermString };                    // Kr_p
    ViewKey dPhaseRelPerm_dPhaseVolFraction = { dPhaseRelPerm_dPhaseVolFractionString }; // dKr_p/dS_p

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
  array3d<real64>  m_phaseRelPerm;
  array4d<real64>  m_dPhaseRelPerm_dPhaseVolFrac;


};


template< typename LEAFCLASS, typename ... ARGS >
void RelativePermeabilityBase::BatchUpdateKernel( arrayView2d<real64 const> const & phaseVolumeFraction,
                                                  ARGS&& ... args)
{
  localIndex const numElem = m_phaseRelPerm.size(0);
  localIndex const numQ = m_phaseRelPerm.size(1);
  localIndex const NP = numFluidPhases();

  arrayView3d<real64> const & phaseRelPerm = m_phaseRelPerm;
  arrayView4d<real64> const & dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac;

  forall_in_range( 0, numElem, GEOSX_LAMBDA ( localIndex const k )
  {
    for( localIndex q=0 ; q<numQ ; ++q )
    {
      LEAFCLASS::StateUpdatePointRelPerm( NP,
                                          phaseVolumeFraction[k],
                                          phaseRelPerm[k][q],
                                          dPhaseRelPerm_dPhaseVolFrac[k][q],
                                          args...);
    }
  });
}

} // namespace constitutive

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
