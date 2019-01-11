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
  * @file SingleFluidBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SINGLEFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SINGLEFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{

class SingleFluidBase : public ConstitutiveBase
{
public:

  SingleFluidBase( std::string const & name, ManagedGroup * const parent );

  virtual ~SingleFluidBase() override;

  // *** ConstitutiveBase interface

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** SingleFluid-specific interface

  /**
   * @brief Perform a single point constitutive update.
   * @param[in] pressure the target pressure value
   * @param[in] k first constitutive index (e.g. elem index)
   * @param[in] q second constitutive index (e.g. quadrature index)
   *
   * @note This function should generally not be called from a kernel, use BatchUpdate instead
   */
  virtual void PointUpdate( real64 const & pressure, localIndex const k, localIndex const q ) = 0;

  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] pressure array containing target pressure values
   */
  virtual void BatchUpdate( arrayView1d<real64 const> const & pressure ) = 0;

  /**
   * @brief Compute constitutive values at a single point.
   * @param[in]  pressure the target pressure value
   * @param[out] density fluid density
   * @param[out] dDensity_dPressure fluid density derivative w.r.t. pressure
   * @param[out] viscosity fluid viscosity
   * @param[out] dViscosity_dPressure fluid viscosity derivatife w.r.t. pressure
   *
   * @note This function should only be called in extremely rare cases, when constitutive state
   * needs to be evaluated at a point where constitutive model does not have storage allocated.
   * It should not be called from kernels since it is virtual.
   */
  virtual void Compute( real64 const & pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const = 0;

  array2d<real64> const & density() const { return m_density; }
  array2d<real64>       & density()       { return m_density; }

  array2d<real64> const & dPressure_dDensity() const { return m_dDensity_dPressure; }
  array2d<real64>       & dPressure_dDensity()       { return m_dDensity_dPressure; }

  array2d<real64> const & viscosity() const { return m_viscosity; }
  array2d<real64>       & viscosity()       { return m_viscosity; }

  array2d<real64> const & dViscosity_dDensity() const { return m_dViscosity_dPressure; }
  array2d<real64>       & dViscosity_dDensity()       { return m_dViscosity_dPressure; }

  // *** Data repository keys

  struct viewKeyStruct
  {
    static constexpr auto densityString      = "density";
    static constexpr auto dDens_dPresString  = "dPressure_dDensity";

    static constexpr auto viscosityString    = "viscosity";
    static constexpr auto dVisc_dPresString  = "dViscosity_dDensity";

    using ViewKey = dataRepository::ViewKey;

    ViewKey density     = { densityString };
    ViewKey dDens_dPres = { dDens_dPresString };

    ViewKey viscosity   = { viscosityString };
    ViewKey dVisc_dPres = { dVisc_dPresString };

  } viewKeysSingleFluidBase;

protected:

  virtual void PostProcessInput() override;

  /**
   * Function to take care of launching the kernel over all constitutive points
   * @tparam POLICY execution policy to use for the launch
   * @tparam LAMBDA type the target function
   * @param lambda the kernel function
   */
  template< typename POLICY=elemPolicy, typename LAMBDA >
  void LaunchKernel( LAMBDA && lambda );

  /**
   * @brief Function to batch process constitutive updates via a kernel launch.
   * @tparam LEAFCLASS The derived class that provides the functions for use in the kernel
   * @tparam POLICY Kernel launch policy (defaults to element policy, but can be chosen by the implementation)
   * @tparam ARGS Parameter pack for arbitrary number of arbitrary types for the function parameter list
   * @param pressure array containing the pressure values,  which is input to the update.
   * @param args arbitrary number of arbitrary types that are passed to the kernel
   *
   * @note This function expects LEAFCLASS to have a public static function Compute with the appropriate signature
   */
  template< typename LEAFCLASS, typename POLICY=elemPolicy, typename ... ARGS >
  void BatchUpdateKernel( arrayView1d<real64 const> const & pressure,
                          ARGS && ... args );

  /**
   * @brief Function to batch process density updates via a kernel launch.
   * @tparam LEAFCLASS The derived class that provides the functions for use in the kernel
   * @tparam POLICY Kernel launch policy (defaults to element policy, but can be chosen by the implementation)
   * @tparam ARGS Parameter pack for arbitrary number of arbitrary types for the function parameter list
   * @param pressure array containing the pressure values,  which is input to the update.
   * @param args arbitrary number of arbitrary types that are passed to the kernel
   *
   * @note This function expects LEAFCLASS to have a public static function Compute with the appropriate signature
   */
  template< typename LEAFCLASS, typename POLICY=elemPolicy, typename ... ARGS >
  void BatchDensityUpdateKernel( arrayView1d<real64 const> const & pressure,
                                 ARGS && ... args );

  /**
   * @brief Function to batch process viscosity updates via a kernel launch.
   * @tparam LEAFCLASS The derived class that provides the functions for use in the kernel
   * @tparam POLICY Kernel launch policy (defaults to element policy, but can be chosen by the implementation)
   * @tparam ARGS Parameter pack for arbitrary number of arbitrary types for the function parameter list
   * @param pressure array containing the pressure values,  which is input to the update.
   * @param args arbitrary number of arbitrary types that are passed to the kernel
   *
   * @note This function expects LEAFCLASS to have a public static function Compute with the appropriate signature
   */
  template< typename LEAFCLASS, typename POLICY=elemPolicy, typename ... ARGS >
  void BatchViscosityUpdateKernel( arrayView1d<real64 const> const & pressure,
                                   ARGS && ... args );

  array2d<real64> m_density;
  array2d<real64> m_dDensity_dPressure;

  array2d<real64> m_viscosity;
  array2d<real64> m_dViscosity_dPressure;

};


template< typename POLICY, typename LAMBDA >
void SingleFluidBase::LaunchKernel( LAMBDA && lambda )
{
  localIndex const numElem = m_density.size(0);
  localIndex const numQuad = m_density.size(1);

  forall_in_range<POLICY>( 0, numElem, GEOSX_LAMBDA ( localIndex const k )
  {
    for (localIndex q = 0; q < numQuad; ++q)
    {
      lambda( k, q );
    }
  } );
}

template< typename LEAFCLASS, typename POLICY, typename ... ARGS >
void SingleFluidBase::BatchUpdateKernel( arrayView1d<real64 const> const & pressure,
                                         ARGS && ... args )
{
  arrayView2d<real64> const & density = m_density;
  arrayView2d<real64> const & dDensity_dPressure = m_dDensity_dPressure;
  arrayView2d<real64> const & viscosity = m_viscosity;
  arrayView2d<real64> const & dViscosity_dPressure = m_dViscosity_dPressure;

  LaunchKernel<POLICY>( GEOSX_LAMBDA ( localIndex const k, localIndex const q )
  {
    LEAFCLASS::Compute( pressure[k],
                        density[k][q],
                        dDensity_dPressure[k][q],
                        viscosity[k][q],
                        dViscosity_dPressure[k][q],
                        args... );
  } );
}

template< typename LEAFCLASS, typename POLICY, typename ... ARGS >
void SingleFluidBase::BatchDensityUpdateKernel( arrayView1d<real64 const> const & pressure,
                                                ARGS && ... args )
{
  arrayView2d<real64> const & density = m_density;
  arrayView2d<real64> const & dDensity_dPressure = m_dDensity_dPressure;

  LaunchKernel<POLICY>( GEOSX_LAMBDA ( localIndex const k, localIndex const q )
  {
    LEAFCLASS::Compute( pressure[k],
                        density[k][q],
                        dDensity_dPressure[k][q],
                        args... );
  } );
}

template< typename LEAFCLASS, typename POLICY, typename ... ARGS >
void SingleFluidBase::BatchViscosityUpdateKernel( arrayView1d<real64 const> const & pressure,
                                                  ARGS && ... args )
{
  arrayView2d<real64> const & viscosity = m_viscosity;
  arrayView2d<real64> const & dViscosity_dPressure = m_dViscosity_dPressure;

  LaunchKernel<POLICY>( GEOSX_LAMBDA ( localIndex const k, localIndex const q )
  {
    LEAFCLASS::Compute( pressure[k],
                        viscosity[k][q],
                        dViscosity_dPressure[k][q],
                        args... );
  } );
}

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SINGLEFLUIDBASE_HPP
