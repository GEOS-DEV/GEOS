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
 * @file SlurryFluidBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class SlurryFluidBase
 * A class to calculate slurry fluid properties
 */

class SlurryFluidBase : public ConstitutiveBase
{
public:

  SlurryFluidBase( std::string const & name, Group * const parent );

  virtual ~SlurryFluidBase() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const override = 0;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  static constexpr localIndex MAX_NUM_COMPONENTS = 32;

  // *** SlurryFluidBase-specific interface

  /**
   * @brief Perform a single point constitutive update.
   * @param[in] pressure, pressure value
   * @param[in] proppantConcentration, proppant concentration value
   * @param[in] Componentconcentration, fluid composition array
   * @param[in] shearRate, shear rate for power-law fluid calculation
   * @param[in] isProppantBoundary, proppant boundary flag
   * @param[in] k, first constitutive index (e.g. elem index)
   * @param[in] q, second constitutive index (e.g. quadrature index)
   *
   * @note This function should generally not be called from a kernel, use BatchUpdate
   */

  virtual void PointUpdate( real64 const & pressure,
                            real64 const & proppantConcentration,
                            arraySlice1d< real64 const > const & Componentconcentration,
                            real64 const & shearRate,
                            integer const & isProppantBoundary,
                            localIndex const k,
                            localIndex const q ) = 0;

  /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] pressure, pressure array containing target pressure values
   * @param[in] proppantConcentration, proppant concentration array containing target pressure values
   * @param[in] componentConcentration, fluid composition array containing target pressure values
   * @param[in] shearRate, shear rate array containing target pressure values
   *
   * @note This function is not implemented
   */
  virtual void BatchUpdate( arrayView1d< real64 const > const & pressure,
                            arrayView1d< real64 const > const & proppantConcentration,
                            arrayView2d< real64 const > const & componentConcentration,
                            arrayView1d< real64 const > const & shearRate ) = 0;

  /**
   * @brief Perform a single point fluid property update.
   * @param[in] pressure, pressure value
   * @param[in] componentconcentration, fluid composition array
   * @param[in] shearRate, shear rate for power-law fluid calculation
   * @param[in] k, first constitutive index (e.g. elem index)
   * @param[in] q, second constitutive index (e.g. quadrature index)
   *
   */
  virtual void PointUpdateFluidProperty( real64 const & pressure,
                                         arraySlice1d< real64 const > const & componentconcentration,
                                         real64 const & shearRate,
                                         localIndex const k,
                                         localIndex const q ) = 0;

  /**
   * @brief Perform a single point fluid component property update.
   * @param[in] pressure, pressure value
   * @param[in] componentconcentration, fluid composition array
   * @param[in] k, first constitutive index (e.g. elem index)
   * @param[in] q, second constitutive index (e.g. quadrature index)
   *
   */
  virtual void PointUpdateComponentDensity( real64 const & pressure,
                                            arraySlice1d< real64 const > const & componentconcentration,
                                            localIndex const k,
                                            localIndex const q ) = 0;

  localIndex numFluidComponents() const;

  string const & componentName( localIndex ic ) const;

  arrayView1d< real64 const > const & nIndex() const { return m_nIndices; }
  arrayView1d< real64 > const & KIndex() const { return m_Ks; }

  arrayView2d< real64 const > const & density() const { return m_density; }
  arrayView2d< real64 > const & density()       { return m_density; }

  arrayView2d< real64 const > const & dDensity_dPressure() const { return m_dDens_dPres; }
  arrayView2d< real64 > const & dDensity_dPressure()       { return m_dDens_dPres; }

  arrayView2d< real64 const > const & dDensity_dProppantConcentration() const { return m_dDens_dProppantConc; }
  arrayView2d< real64 > const & dDensity_dProppantConcentration()       { return m_dDens_dProppantConc; }

  arrayView3d< real64 const > const & dDensity_dComponentConcentration() const { return m_dDens_dCompConc; }
  arrayView3d< real64 > const & dDensity_dComponentConcentration()       { return m_dDens_dCompConc; }

  arrayView2d< real64 const > const & fluidDensity() const { return m_fluidDensity; }
  arrayView2d< real64 > const & fluidDensity()       { return m_fluidDensity; }

  arrayView2d< real64 const > const & dFluidDensity_dPressure() const { return m_dFluidDens_dPres; }
  arrayView2d< real64 > const & dFluidDensity_dPressure()       { return m_dFluidDens_dPres; }

  arrayView3d< real64 const > const & dFluidDensity_dComponentConcentration() const { return m_dFluidDens_dCompConc; }
  arrayView3d< real64 > const & dFluidDensity_dComponentConcentration()       { return m_dFluidDens_dCompConc; }


  arrayView2d< real64 const > const & fluidViscosity() const { return m_fluidViscosity; }
  arrayView2d< real64 > const & fluidViscosity()       { return m_fluidViscosity; }

  arrayView2d< real64 const > const & dFluidViscosity_dPressure() const { return m_dFluidVisc_dPres; }
  arrayView2d< real64 > const & dFluidViscosity_dPressure()       { return m_dFluidVisc_dPres; }

  arrayView3d< real64 const > const & dFluidViscosity_dComponentConcentration() const { return m_dFluidVisc_dCompConc; }
  arrayView3d< real64 > const & dFluidViscosity_dComponentConcentration()       { return m_dFluidVisc_dCompConc; }


  arrayView2d< real64 const > const & viscosity() const { return m_viscosity; }
  arrayView2d< real64 > const & viscosity()       { return m_viscosity; }

  arrayView2d< real64 const > const & dViscosity_dPressure() const { return m_dVisc_dPres; }
  arrayView2d< real64 > const & dViscosity_dPressure()       { return m_dVisc_dPres; }

  arrayView2d< real64 const > const & dViscosity_dProppantConcentration() const { return m_dVisc_dProppantConc; }
  arrayView2d< real64 > const & dViscosity_dProppantConcentration()       { return m_dVisc_dProppantConc; }

  arrayView3d< real64 const > const & dViscosity_dComponentConcentration() const { return m_dVisc_dCompConc; }
  arrayView3d< real64 > const & dViscosity_dComponentConcentration()       { return m_dVisc_dCompConc; }

  bool isNewtonianFluid() const { return m_isNewtonianFluid; }

  // *** Data repository keys

  struct viewKeyStruct
  {

    static constexpr auto componentNamesString       = "componentNames";
    static constexpr auto defaultDensityString      = "defaultDensity";
    static constexpr auto defaultCompressibilityString      = "defaultCompressibility";
    static constexpr auto defaultViscosityString      = "defaultViscosity";

    static constexpr auto densityString      = "density";
    static constexpr auto dDens_dPresString  = "dDens_dPres";
    static constexpr auto dDens_dProppantConcString  = "dDens_dProppantConc";
    static constexpr auto dDens_dCompConcString  = "dDens_dCompConc";

    static constexpr auto componentDensityString      = "componentDensity";
    static constexpr auto dCompDens_dPresString  = "dCompDens_dPres";
    static constexpr auto dCompDens_dCompConcString  = "dCompDens_dCompConc";

    static constexpr auto fluidDensityString      = "FluidDensity";
    static constexpr auto dFluidDens_dPresString  = "dFluidDens_dPres";
    static constexpr auto dFluidDens_dCompConcString  = "dFluidDens_dCompConc";

    static constexpr auto fluidViscosityString      = "FluidViscosity";
    static constexpr auto dFluidVisc_dPresString  = "dFluidVisc_dPres";
    static constexpr auto dFluidVisc_dCompConcString  = "dFluidVisc_dCompConc";

    static constexpr auto viscosityString    = "viscosity";
    static constexpr auto dVisc_dPresString  = "dVisc_dPres";
    static constexpr auto dVisc_dProppantConcString  = "dVisc_dProppantConc";
    static constexpr auto dVisc_dCompConcString  = "dVisc_dCompConc";
    static constexpr auto flowBehaviorIndexString   = "flowBehaviorIndex";

    static constexpr auto flowConsistencyIndexString   = "flowConsistencyIndex";
  } viewKeysSlurryFluidBase;

protected:

  virtual void PostProcessInput() override;

  string_array m_componentNames;

  array1d< real64 > m_defaultDensity;
  array1d< real64 > m_defaultCompressibility;
  array1d< real64 > m_defaultViscosity;

  array1d< real64 > m_nIndices;
  array1d< real64 > m_Ks;

  array2d< real64 > m_density;
  array2d< real64 > m_dDens_dPres;
  array2d< real64 > m_dDens_dProppantConc;
  array3d< real64 > m_dDens_dCompConc;

  array3d< real64 > m_componentDensity;
  array3d< real64 > m_dCompDens_dPres;
  array4d< real64 > m_dCompDens_dCompConc;

  array2d< real64 > m_fluidDensity;
  array2d< real64 > m_dFluidDens_dPres;
  array3d< real64 > m_dFluidDens_dCompConc;

  array2d< real64 > m_fluidViscosity;
  array2d< real64 > m_dFluidVisc_dPres;
  array3d< real64 > m_dFluidVisc_dCompConc;

  array2d< real64 > m_viscosity;
  array2d< real64 > m_dVisc_dPres;
  array2d< real64 > m_dVisc_dProppantConc;
  array3d< real64 > m_dVisc_dCompConc;

  bool m_isNewtonianFluid;

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_SLURRYFLUIDBASE_HPP
