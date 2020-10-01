/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProppantSlurryFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PROPPANTSLURRYFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PROPPANTSLURRYFLUID_HPP_

#include "constitutive/fluid/SlurryFluidBase.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @brief Update class for the proppant slurry model suitable for lambda capture.
 */
class ProppantSlurryFluidUpdate final : public SlurryFluidBaseUpdate
{
public:

  /**
   * @brief
   * @param compressibility
   * @param referenceProppantDensity
   * @param referencePressure
   * @param referenceDensity
   * @param referenceViscosity
   * @param maxProppantConcentration
   * @param defaultDensity
   * @param defaultCompressibility
   * @param defaultViscosity
   * @param nIndices
   * @param Ks
   * @param isNewtonianFluid
   * @param density
   * @param dDens_dPres
   * @param dDens_dProppantConc
   * @param dDens_dCompConc
   * @param componentDensity
   * @param dCompDens_dPres
   * @param dCompDens_dCompConc
   * @param fluidDensity
   * @param dFluidDens_dPres
   * @param dFluidDens_dCompConc
   * @param fluidViscosity
   * @param dFluidVisc_dPres
   * @param dFluidVisc_dCompConc
   * @param viscosity
   * @param dVisc_dPres
   * @param dVisc_dProppantConc
   * @param dVisc_dCompConc
   */
  ProppantSlurryFluidUpdate( real64 const compressibility,
                             real64 const referenceProppantDensity,
                             real64 const referencePressure,
                             real64 const referenceDensity,
                             real64 const referenceViscosity,
                             real64 const maxProppantConcentration,
                             arrayView1d< real64 const > const & defaultDensity,
                             arrayView1d< real64 const > const & defaultCompressibility,
                             arrayView1d< real64 const > const & defaultViscosity,
                             arrayView1d< real64 const > const & nIndices,
                             arrayView1d< real64 const > const & ks,
                             bool const isNewtonianFluid,
                             arrayView2d< real64 > const & density,
                             arrayView2d< real64 > const & dDensDPres,
                             arrayView2d< real64 > const & dDensDProppantConc,
                             arrayView3d< real64 > const & dDensDCompConc,
                             arrayView3d< real64 > const & componentDensity,
                             arrayView3d< real64 > const & dCompDensDPres,
                             arrayView4d< real64 > const & dCompDensDCompConc,
                             arrayView2d< real64 > const & fluidDensity,
                             arrayView2d< real64 > const & dFluidDensDPres,
                             arrayView3d< real64 > const & dFluidDensDCompConc,
                             arrayView2d< real64 > const & fluidViscosity,
                             arrayView2d< real64 > const & dFluidViscDPres,
                             arrayView3d< real64 > const & dFluidViscDCompConc,
                             arrayView2d< real64 > const & viscosity,
                             arrayView2d< real64 > const & dViscDPres,
                             arrayView2d< real64 > const & dViscDProppantConc,
                             arrayView3d< real64 > const & dViscDCompConc )
    : SlurryFluidBaseUpdate( defaultDensity,
                             defaultCompressibility,
                             defaultViscosity,
                             nIndices,
                             ks,
                             isNewtonianFluid,
                             density,
                             dDensDPres,
                             dDensDProppantConc,
                             dDensDCompConc,
                             componentDensity,
                             dCompDensDPres,
                             dCompDensDCompConc,
                             fluidDensity,
                             dFluidDensDPres,
                             dFluidDensDCompConc,
                             fluidViscosity,
                             dFluidViscDPres,
                             dFluidViscDCompConc,
                             viscosity,
                             dViscDPres,
                             dViscDProppantConc,
                             dViscDCompConc ),
    m_compressibility( compressibility ),
    m_referenceProppantDensity( referenceProppantDensity ),
    m_referencePressure( referencePressure ),
    m_referenceDensity( referenceDensity ),
    m_referenceViscosity( referenceViscosity ),
    m_maxProppantConcentration( maxProppantConcentration )
  {}

  /// Default copy constructor
  ProppantSlurryFluidUpdate( ProppantSlurryFluidUpdate const & ) = default;

  /// Default move constructor
  ProppantSlurryFluidUpdate( ProppantSlurryFluidUpdate && ) = default;

  /// Deleted copy assignment operator
  ProppantSlurryFluidUpdate & operator=( ProppantSlurryFluidUpdate const & ) = delete;

  /// Deleted move assignment operator
  ProppantSlurryFluidUpdate & operator=( ProppantSlurryFluidUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void Update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const proppantConcentration,
                       arraySlice1d< real64 const > const & componentConcentration,
                       real64 const shearRate,
                       integer const isProppantBoundary ) const override
  {
    UpdateFluidProperty( k,
                         q,
                         pressure,
                         componentConcentration,
                         shearRate );

    Compute( proppantConcentration,
             m_fluidDensity[k][q],
             m_dFluidDens_dPres[k][q],
             m_dFluidDens_dCompConc[k][q],
             m_fluidViscosity[k][q],
             m_dFluidVisc_dPres[k][q],
             m_dFluidVisc_dCompConc[k][q],
             isProppantBoundary,
             m_density[k][q],
             m_dDens_dPres[k][q],
             m_dDens_dProppantConc[k][q],
             m_dDens_dCompConc[k][q],
             m_viscosity[k][q],
             m_dVisc_dPres[k][q],
             m_dVisc_dProppantConc[k][q],
             m_dVisc_dCompConc[k][q] );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void UpdateFluidProperty( localIndex const k,
                                    localIndex const q,
                                    real64 const pressure,
                                    arraySlice1d< real64 const > const & componentConcentration,
                                    real64 const shearRate ) const override
  {
    GEOSX_UNUSED_VAR( shearRate )
    ComputeFluidDensity( pressure,
                         componentConcentration,
                         m_componentDensity[k][q],
                         m_dCompDens_dPres[k][q],
                         m_dCompDens_dCompConc[k][q],
                         m_fluidDensity[k][q],
                         m_dFluidDens_dPres[k][q],
                         m_dFluidDens_dCompConc[k][q] );

    ComputeFluidViscosity( m_componentDensity[k][q],
                           m_dCompDens_dPres[k][q],
                           m_dCompDens_dCompConc[k][q],
                           m_fluidDensity[k][q],
                           m_dFluidDens_dPres[k][q],
                           m_dFluidDens_dCompConc[k][q],
                           m_fluidViscosity[k][q],
                           m_dFluidVisc_dPres[k][q],
                           m_dFluidVisc_dCompConc[k][q] );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void UpdateComponentDensity( localIndex const k,
                                       localIndex const q,
                                       real64 const pressure,
                                       arraySlice1d< real64 const > const & componentConcentration ) const override
  {
    ComputeComponentDensity( pressure,
                             componentConcentration,
                             m_componentDensity[k][q],
                             m_dCompDens_dPres[k][q],
                             m_dCompDens_dCompConc[k][q] );
  }

private:

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void ComputeFluidDensity( real64 const & pressure,
                            arraySlice1d< real64 const > const & componentConcentration,
                            arraySlice1d< real64 > const & componentDensity,
                            arraySlice1d< real64 > const & dComponentDensityDPressure,
                            arraySlice2d< real64 > const & dComponentDensityDComponentConcentration,
                            real64 & fluidDensity,
                            real64 & dFluidDensityDPressure,
                            arraySlice1d< real64 > const & dFluidDensityDComponentConcentration ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void ComputeComponentDensity( real64 const pressure,
                                arraySlice1d< real64 const > const & componentConcentration,
                                arraySlice1d< real64 > const & componentDensity,
                                arraySlice1d< real64 > const & dComponentDensityDPressure,
                                arraySlice2d< real64 > const & dComponentDensityDComponentConcentration ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void ComputeFluidViscosity( arraySlice1d< real64 const > const & componentDensity,
                              arraySlice1d< real64 const > const & dComponentDensityDPressure,
                              arraySlice2d< real64 const > const & dComponentDensityDComponentConcentration,
                              real64 const fluidDensity,
                              real64 const dFluidDensityDPressure,
                              arraySlice1d< real64 const > const & dFluidDensityDComponentConcentration,
                              real64 & fluidViscosity,
                              real64 & dFluidViscosityDPressure,
                              arraySlice1d< real64 > const & dFluidViscosityDComponentConcentration ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void Compute( real64 const & proppantConcentration,
                real64 const & fluidDensity,
                real64 const & dFluidDensityDPressure,
                arraySlice1d< real64 const > const & dFluidDensityDComponentConcentration,
                real64 const & fluidViscosity,
                real64 const & dFluidViscosityDPressure,
                arraySlice1d< real64 const > const & dFluidViscosityDComponentConcentration,
                integer const & isProppantBoundary,
                real64 & density,
                real64 & dDensityDPressure,
                real64 & dDensityDProppantConcentration,
                arraySlice1d< real64 > const & dDensityDComponentConcentration,
                real64 & viscosity,
                real64 & dViscosityDPressure,
                real64 & dViscosityDProppantConcentration,
                arraySlice1d< real64 > const & dViscosityDComponentConcentration ) const;

  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantConcentration;

};

class ProppantSlurryFluid : public SlurryFluidBase
{
public:

  ProppantSlurryFluid( std::string const & name, Group * const parent );

  virtual ~ProppantSlurryFluid() override;

  // *** ConstitutiveBase interface
  static std::string CatalogName() { return "ProppantSlurryFluid"; }

  virtual string getCatalogName() const override { return CatalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  using KernelWrapper = ProppantSlurryFluidUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  // *** Data repository keys

  struct viewKeyStruct : public SlurryFluidBase::viewKeyStruct
  {
    static constexpr auto compressibilityString    = "compressibility";
    static constexpr auto referencePressureString  = "referencePressure";
    static constexpr auto referenceDensityString   = "referenceDensity";
    static constexpr auto referenceProppantDensityString   = "referenceProppantDensity";
    static constexpr auto maxProppantConcentrationString   = "maxProppantConcentration";
    static constexpr auto referenceViscosityString = "referenceViscosity";

  } viewKeysProppantSlurryFluid;

protected:

  virtual void PostProcessInput() override;

private:

  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantConcentration;

};

GEOSX_HOST_DEVICE
void
ProppantSlurryFluidUpdate::
  ComputeFluidDensity( real64 const & pressure,
                       arraySlice1d< real64 const > const & componentConcentration,
                       arraySlice1d< real64 > const & componentDensity,
                       arraySlice1d< real64 > const & dComponentDensityDPressure,
                       arraySlice2d< real64 > const & dComponentDensityDComponentConcentration,
                       real64 & fluidDensity,
                       real64 & dFluidDensityDPressure,
                       arraySlice1d< real64 > const & dFluidDensityDComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  real64 const baseFluidDensity = m_referenceDensity * exp( m_compressibility * (pressure - m_referencePressure));
  real64 const dBaseFluidDensity_dPressure = m_compressibility * baseFluidDensity;

  fluidDensity = baseFluidDensity;
  dFluidDensityDPressure = dBaseFluidDensity_dPressure;
  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidDensityDComponentConcentration[c] = 0.0;
  }

  for( localIndex c = 0; c < NC; ++c )
  {
    real64 const density = m_defaultDensity[c] * exp( m_defaultCompressibility[c] * (pressure - m_referencePressure));

    componentDensity[c] = componentConcentration[c] * density;
    dComponentDensityDPressure[c] = m_defaultCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
    {
      dComponentDensityDComponentConcentration[c][i] = 0.0;
    }
    dComponentDensityDComponentConcentration[c][c] = density;

    fluidDensity += componentDensity[c] - componentConcentration[c] * baseFluidDensity;
    dFluidDensityDPressure += dComponentDensityDPressure[c] - componentConcentration[c] * dBaseFluidDensity_dPressure;
    dFluidDensityDComponentConcentration[c] += dComponentDensityDComponentConcentration[c][c] - baseFluidDensity;
  }

  // TODO: why? - Sergey
  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidDensityDComponentConcentration[c] = 0.0;
  }
}

GEOSX_HOST_DEVICE
void
ProppantSlurryFluidUpdate::
  ComputeComponentDensity( real64 const pressure,
                           arraySlice1d< real64 const > const & componentConcentration,
                           arraySlice1d< real64 > const & componentDensity,
                           arraySlice1d< real64 > const & dComponentDensityDPressure,
                           arraySlice2d< real64 > const & dComponentDensityDComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  for( localIndex c = 0; c < NC; ++c )
  {
    real64 const density = m_defaultDensity[c] * exp( m_defaultCompressibility[c] * (pressure - m_referencePressure));

    componentDensity[c] = componentConcentration[c] * density;
    dComponentDensityDPressure[c] = m_defaultCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
    {
      dComponentDensityDComponentConcentration[c][i] = 0.0;
    }
    dComponentDensityDComponentConcentration[c][c] = density;
  }
}

GEOSX_HOST_DEVICE
void
ProppantSlurryFluidUpdate::
  ComputeFluidViscosity( arraySlice1d< real64 const > const & componentDensity,
                         arraySlice1d< real64 const > const & dComponentDensityDPressure,
                         arraySlice2d< real64 const > const & GEOSX_UNUSED_PARAM( dComponentDensity_dComponentConcentration ),
                         real64 const fluidDensity,
                         real64 const dFluidDensityDPressure,
                         arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
                         real64 & fluidViscosity,
                         real64 & dFluidViscosityDPressure,
                         arraySlice1d< real64 > const & dFluidViscosityDComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  fluidViscosity = m_referenceViscosity;
  dFluidViscosityDPressure = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidViscosityDComponentConcentration[c] = 0.0;
  }

  for( localIndex c1 = 0; c1 < NC; ++c1 )
  {
    fluidViscosity += componentDensity[c1] / fluidDensity * (m_defaultViscosity[c1] - m_referenceViscosity);

    dFluidViscosityDPressure +=
      (dComponentDensityDPressure[c1] / fluidDensity - componentDensity[c1] / fluidDensity / fluidDensity * dFluidDensityDPressure) *
      (m_defaultViscosity[c1] - m_referenceViscosity);
  }
}

void
ProppantSlurryFluidUpdate::
  Compute( real64 const & proppantConcentration,
           real64 const & fluidDensity,
           real64 const & dFluidDensityDPressure,
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
           real64 const & fluidViscosity,
           real64 const & dFluidViscosityDPressure,
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidViscosity_dComponentConcentration ),
           integer const & isProppantBoundary,
           real64 & density,
           real64 & dDensityDPressure,
           real64 & dDensityDProppantConcentration,
           arraySlice1d< real64 > const & dDensityDComponentConcentration,
           real64 & viscosity,
           real64 & dViscosityDPressure,
           real64 & dViscosityDProppantConcentration,
           arraySlice1d< real64 > const & dViscosityDComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  real64 effectiveConcentration = proppantConcentration;
  if( effectiveConcentration > 0.95 * m_maxProppantConcentration || isProppantBoundary == 1 )
  {
    effectiveConcentration = 0.0;
  }

  density = (1.0 - effectiveConcentration) * fluidDensity + effectiveConcentration * m_referenceProppantDensity;
  dDensityDPressure = (1.0 - effectiveConcentration) * dFluidDensityDPressure;

  dDensityDProppantConcentration = 0.0;
  for( localIndex c = 0; c < NC; ++c )
  {
    dDensityDComponentConcentration[c] = 0.0;
  }

  real64 const coef = pow( 1.0 + 1.25 *  effectiveConcentration / (1.0 - effectiveConcentration / m_maxProppantConcentration), 2.0 );
  viscosity = fluidViscosity * coef;
  dViscosityDPressure = dFluidViscosityDPressure * coef;

  dViscosityDProppantConcentration = 0.0;
  for( localIndex c = 0; c < NC; ++c )
  {
    dViscosityDComponentConcentration[c] = 0.0;
  }
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_PROPPANTSLURRYFLUID_HPP_ */
