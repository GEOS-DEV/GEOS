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
 *  @file DruckerPragerExtended.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP_

#include "ElastoPlasticUpdates.hpp"
#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class DruckerPragerExtendedUpdates
 *
 * Class to provide material updates that may be called from a kernel function.
 */
class DruckerPragerExtendedUpdates : public ElastoPlasticUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] initialFriction The ArrayView holding the initial friction data for each element.
   * @param[in] residualFriction The ArrayView holding the residual friction data for each element.
   * @param[in] dilationRatio The ArrayView holding the ratio between dilation and friction data for each element.
   * @param[in] pressureIntercept The ArrayView holding the pressure intercept for each element.
   * @param[in] hardening The ArrayView holding the hardening data for each element.
   * @param[in] state The ArrayView holding the state data for each element.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  DruckerPragerExtendedUpdates( arrayView1d< real64 const > const & initialFriction,
                                arrayView1d< real64 const > const & residualFriction,
                                arrayView1d< real64 const > const & dilationRatio,
                                arrayView1d< real64 const > const & pressureIntercept,
                                arrayView1d< real64 const > const & hardening,
                                arrayView2d< real64 > const & state,
                                arrayView1d< real64 const > const & bulkModulus,
                                arrayView1d< real64 const > const & shearModulus,
                                arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    ElastoPlasticUpdates( state, bulkModulus, shearModulus, newStress, oldStress ),
    m_initialFriction( initialFriction ),
    m_residualFriction( residualFriction ),
    m_dilationRatio( dilationRatio ),
    m_pressureIntercept( pressureIntercept ),
    m_hardening( hardening ),
    m_state( state )
  {}

  /// Default copy constructor
  DruckerPragerExtendedUpdates( DruckerPragerExtendedUpdates const & ) = default;

  /// Default move constructor
  DruckerPragerExtendedUpdates( DruckerPragerExtendedUpdates && ) = default;

  /// Deleted default constructor
  DruckerPragerExtendedUpdates() = delete;

  /// Deleted copy assignment operator
  DruckerPragerExtendedUpdates & operator=( DruckerPragerExtendedUpdates const & ) = delete;

  /// Deleted move assignment operator
  DruckerPragerExtendedUpdates & operator=( DruckerPragerExtendedUpdates && ) =  delete;

private:

  /// A reference to the ArrayView holding the initial friction coefficient for each element.
  arrayView1d< real64 const > const m_initialFriction;

  /// A reference to the ArrayView holding the residual friction coefficient for each element.
  arrayView1d< real64 const > const m_residualFriction;

  /// A reference to the ArrayView holding the dilation ratio for each element.
  arrayView1d< real64 const > const m_dilationRatio;

  /// A reference to the ArrayView holding the pressure intercept for each element.
  arrayView1d< real64 const > const m_pressureIntercept;

  /// A reference to the ArrayView holding the hardening parameter for each element.
  arrayView1d< real64 const > const m_hardening;

  /// A reference to the ArrayView holding the state variable for each integration point
  arrayView2d< real64 > const m_state;

  virtual real64 yield( localIndex const k,
                        localIndex const GEOSX_UNUSED_PARAM( q ),
                        real64 const invP,
                        real64 const invQ,
                        real64 const friction ) const override final
  {
    return  invQ + friction * ( invP - m_pressureIntercept[k] );
  }

  virtual void yieldDerivatives( localIndex const k,
                                 localIndex const GEOSX_UNUSED_PARAM( q ),
                                 real64 const invP,
                                 real64 const GEOSX_UNUSED_PARAM( invQ ),
                                 real64 const friction,
                                 real64 (& dF)[3] ) const override final
  {

    // The yield function is: invQ + friction * ( invP - m_pressureIntercept[k] )

    real64 dF_dP = friction;
    real64 dF_dQ = 1.0;
    real64 dF_dFriction = invP - m_pressureIntercept[k];

    dF[0] = dF_dP;
    dF[1] = dF_dQ;
    dF[2] = dF_dFriction;
  }

  virtual void potentialDerivatives( localIndex const k,
                                     localIndex const GEOSX_UNUSED_PARAM( q ),
                                     real64 const GEOSX_UNUSED_PARAM( invP ),
                                     real64 const GEOSX_UNUSED_PARAM( invQ ),
                                     real64 const friction,
                                     real64 (& dG)[8] ) const override final
  { 

    // The plastic potential function is: G = invQ + invP * m_dilationRatio[k] * friction

    real64 dG_dP = m_dilationRatio[k] * friction;
    real64 dG_dQ = 1.0;

    real64 dG_dP_dP = 0.0;
    real64 dG_dP_dQ = 0.0;
    real64 dG_dP_dH = m_dilationRatio[k];

    real64 dG_dQ_dP = 0.0;
    real64 dG_dQ_dQ = 0.0;
    real64 dG_dQ_dH = 0.0;

    dG[0] = dG_dP;
    dG[1] = dG_dQ;

    dG[2] = dG_dP_dP;
    dG[3] = dG_dP_dQ;
    dG[4] = dG_dP_dH;

    dG[5] = dG_dQ_dP;
    dG[6] = dG_dQ_dQ;
    dG[7] = dG_dQ_dH;
  }

  virtual real64 hardening( localIndex const k,
                            localIndex const GEOSX_UNUSED_PARAM( q ),
                            real64 const state ) const override final
  {
    if( state<1e-9 )
    {
      return m_initialFriction[k];
    }
    else
    {
      return m_initialFriction[k] + ( m_residualFriction[k] - m_initialFriction[k] ) * state / ( m_hardening[k] + state );
    }
  }

  virtual real64 hardeningDerivatives( localIndex const k,
                                       localIndex const GEOSX_UNUSED_PARAM( q ),
                                       real64 const state ) const override final
  {
    if( state<1e-9 )
    {
      return 0.0;
    }
    else
    {
      return ( m_residualFriction[k] - m_initialFriction[k] ) * m_hardening[k] / ( m_hardening[k] + state ) / ( m_hardening[k] + state );
    }
  }

  virtual real64 tmpState( localIndex const k,
                           localIndex const q,
                           real64 const GEOSX_UNUSED_PARAM( invP ),
                           real64 const GEOSX_UNUSED_PARAM( invQ ),
                           real64 const plasticMultiplier ) const override final
  {

    // The temporal state variable is updated inside the Newton loops by: state += plasticMultiplier
    // starting from the previous saved state

    return m_state[k][q] + plasticMultiplier;
  }

  virtual real64 stateDerivatives( localIndex const GEOSX_UNUSED_PARAM( k ),
                                   localIndex const GEOSX_UNUSED_PARAM( q ),
                                   real64 const GEOSX_UNUSED_PARAM( invP ),
                                   real64 const GEOSX_UNUSED_PARAM( invQ ),
                                   real64 const GEOSX_UNUSED_PARAM( plasticMultiplier ) ) const override final
  {
    return 1.0;
  }
};

/**
 * @class DruckerPragerExtended
 *
 * Extended Drucker-Prager material model.
 */
class DruckerPragerExtended : public ElasticIsotropic
{
public:

  /// @typedef Alias for DruckerPragerExtendedUpdates
  using KernelWrapper = DruckerPragerExtendedUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DruckerPragerExtended( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DruckerPragerExtended() override;


  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "DruckerPragerExtended";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default initial friction angle
    static constexpr auto defaultInitialFrictionAngleString = "defaultInitialFrictionAngle";

    /// string/key for default initial friction angle
    static constexpr auto defaultResidualFrictionAngleString = "defaultResidualFrictionAngle";

    /// string/key for default dilation angle
    static constexpr auto defaultDilationRatioString = "defaultDilationRatio";

    /// string/key for default hardening rate
    static constexpr auto defaultHardeningString = "defaultHardening";

    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";

    /// string/key for initial friction angle
    static constexpr auto initialFrictionString  = "initialFriction";

    /// string/key for final friction angle
    static constexpr auto residualFrictionString  = "residualFriction";

    /// string/key for dilation angle
    static constexpr auto dilationRatioString  = "dilationRatio";

    /// string/key for pressure intercept
    static constexpr auto pressureInterceptString  = "pressureIntercept";

    /// string/key for cohesion
    static constexpr auto hardeningString  = "hardening";

    /// string/key for state variable
    static constexpr auto stateString  = "stateVariable";
  };

  /**
   * @brief Create a instantiation of the DruckerPragerExtendedUpdate class that refers to the data in this.
   * @return An instantiation of DruckerPragerExtendedUpdate.
   */
  DruckerPragerExtendedUpdates createKernelUpdates() const
  {
    return DruckerPragerExtendedUpdates( m_initialFriction,
                                         m_residualFriction,
                                         m_dilationRatio,
                                         m_pressureIntercept,
                                         m_hardening,
                                         m_state,
                                         m_bulkModulus,
                                         m_shearModulus,
                                         m_newStress,
                                         m_oldStress );
  }

protected:
  virtual void postProcessInput() override;

  /// Material parameter: The default value of the initial yield surface slope
  real64 m_defaultInitialFrictionAngle;

  /// Material parameter: The default value of the final yield surface slope
  real64 m_defaultResidualFrictionAngle;

  /// Material parameter: The default value of the plastic potential slope ratio
  real64 m_defaultDilationRatio;

  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;

  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardening;

  /// Material parameter: The initial yield surface slope param for each element
  array1d< real64 > m_initialFriction;

  /// Material parameter: The final yield surface slope param for each element
  array1d< real64 > m_residualFriction;

  /// Material parameter: The plastic potential slope param (a ratio w.r.t. current yield surface)
  array1d< real64 > m_dilationRatio;

  /// Material parameter: The pressure intercept (location of cone vertex) for each element
  array1d< real64 > m_pressureIntercept;

  /// Material parameter: The hyperbolic hardening parameter for each element
  array1d< real64 > m_hardening;

  /// State variable: The current equivalent plastic shear strain for each quadrature point
  array2d< real64 > m_state;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP_ */

