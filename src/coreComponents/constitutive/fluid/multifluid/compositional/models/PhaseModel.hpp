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

/**
 * @file PhaseModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHASEMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHASEMODEL_HPP_

#include "NullModel.hpp"
#include "ComponentProperties.hpp"
#include "EquationOfState.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief Struct storing the submodels describing the fluid phase behavior.
 * @tparam DENSITY Class describing the density model
 * @tparam VISCOSITY Class describing the viscosity model
 * @tparam ENTHALPY Class describing the enthalpy model
 */
template< typename DENSITY, typename VISCOSITY, typename ENTHALPY >
struct PhaseModel
{
  using Density = DENSITY;
  using Viscosity = VISCOSITY;
  using Enthalpy = ENTHALPY;

  /**
   * @brief Constructor of the struct, in charge of the instantiation of the submodels
   * @param[in] phaseModelName name of the phase model, used only in the instantiation of the submodels
   * @param[in] phaseIndex the index of the phase
   * @param[in] componentProperties Compositional fluid model component parameters
   * @param[in] equationOfState Equation of state parameters
   */
  PhaseModel( string const & phaseModelName,
              ComponentProperties const & componentProperties,
              EquationOfState const & equationOfState,
              integer const phaseIndex ):
    density( phaseModelName + "_" + Density::catalogName(),
             componentProperties,
             equationOfState,
             phaseIndex ),
    viscosity( phaseModelName + "_" + Viscosity::catalogName(),
               componentProperties,
               equationOfState,
               phaseIndex ),
    enthalpy( phaseModelName + "_" + Enthalpy::catalogName(),
              componentProperties,
              equationOfState,
              phaseIndex )
  {}

  /// The phase density model
  Density density;

  /// The phase viscosity model
  Viscosity viscosity;

  /// The phase enthalpy model (can be NoOp for non-thermal models)
  Enthalpy enthalpy;

  /**
   * @brief Struct storing the submodels wrappers used for in-kernel fluid updates
   */
  struct KernelWrapper
  {

    /**
     * @brief Constructor for the kernel wrapper
     * @param[in] dens the density model
     * @param[in] visc the viscosity model
     * @param[in] enth the enthalpy model
     */
    KernelWrapper( Density const & dens,
                   Viscosity const & visc,
                   Enthalpy const & enth ):
      density( dens.createKernelWrapper() ),
      viscosity( visc.createKernelWrapper() ),
      enthalpy( enth.createKernelWrapper() )
    {}

    /**
     * @brief Moves the submodel data to the GPU if the phaseModel is stored in ArrayView
     * @param[in] space the target memory space
     * @param[in] touch the flag to register touch or not
     */
    void move( LvArray::MemorySpace const space, bool const touch )
    {
      density.move( space, touch );
      viscosity.move( space, touch );
      enthalpy.move( space, touch );
    }

    /// Kernel wrapper for density updates
    typename Density::KernelWrapper density;

    /// Kernel wrapper for viscosity updates
    typename Viscosity::KernelWrapper viscosity;

    /// Kernel wrapper for enthalpy updates
    typename Enthalpy::KernelWrapper enthalpy;

  };

  /**
   * @brief Function to create and return a KernelWrapper
   * @return the KernelWrapper object
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( density,
                          viscosity,
                          enthalpy );
  }
};

// A no-op phase model
using NullPhaseModel = PhaseModel< NullModel, NullModel, NullModel >;

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHASEMODEL_HPP_
