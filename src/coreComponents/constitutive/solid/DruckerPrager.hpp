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
 *  @file DruckerPrager.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP_

#include "ElasticIsotropic.hpp"
#include "ElastoPlasticUpdates.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class DruckerPragerUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class DruckerPragerUpdates : public ElastoPlasticUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] friction The ArrayView holding the friction data for each element.
   * @param[in] dilation The ArrayView holding the dilation data for each element.
   * @param[in] hardening The ArrayView holding the hardening data for each element.
   * @param[in] cohesion The ArrayView holding the cohesion data for each element.
   * @param[in] state The ArrayView holding the state data for each element.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  DruckerPragerUpdates( arrayView1d< real64 const > const & friction,
                        arrayView1d< real64 const > const & dilation,
                        arrayView1d< real64 const > const & hardening,
                        arrayView2d< real64 > const & cohesion,
                        arrayView2d< real64 > const & state,
                        arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress ):// TODO tmp stress[6] can be considered 
                                                                                     // in the Elasto-Plastic Newton loops
                                                                                     // to avoid holding both new and old stress 
                                                                                     // on the system
    ElastoPlasticUpdates( bulkModulus, shearModulus, newStress, oldStress ),
    m_friction( friction ),
    m_dilation( dilation ),
    m_hardening( hardening ),
    m_cohesion( cohesion ),
    m_state( state )
  {}

  /// Default copy constructor
  DruckerPragerUpdates( DruckerPragerUpdates const & ) = default;

  /// Default move constructor
  DruckerPragerUpdates( DruckerPragerUpdates && ) = default;

  /// Deleted default constructor
  DruckerPragerUpdates() = delete;

  /// Deleted copy assignment operator
  DruckerPragerUpdates & operator=( DruckerPragerUpdates const & ) = delete;

  /// Deleted move assignment operator
  DruckerPragerUpdates & operator=( DruckerPragerUpdates && ) =  delete;

private:

  virtual real64 yield( localIndex const k,
                        localIndex const GEOSX_UNUSED_PARAM( q ),
                        real64 const invP,
                        real64 const invQ,
                        real64 const cohesion ) const override final
  {
    return  invQ + m_friction[k] * invP - cohesion;
  }

  virtual void yieldDerivatives( localIndex const k,
                                 localIndex const GEOSX_UNUSED_PARAM( q ),
                                 real64 const GEOSX_UNUSED_PARAM( invP ),
                                 real64 const GEOSX_UNUSED_PARAM( invQ ),
                                 real64 const GEOSX_UNUSED_PARAM( cohesion ),
                                 real64 (& dF)[3] ) const override final
  {

    // The yield function is: F = Q + friction * P - cohesion

    real64 dF_dP = m_friction[k];
    real64 dF_dQ = 1.0;
    real64 dF_dCohesion = -1.0;

    dF[0] = dF_dP;
    dF[1] = dF_dQ;
    dF[2] = dF_dCohesion;
  }

  // Derivatives of the plastic potential function to the stress invariants and the hardening parameter.

  virtual void potentialDerivatives( localIndex const k,
                                     localIndex const GEOSX_UNUSED_PARAM( q ),
                                     real64 const GEOSX_UNUSED_PARAM( invP ),
                                     real64 const GEOSX_UNUSED_PARAM( invQ ),
                                     real64 const GEOSX_UNUSED_PARAM( cohesion ),
                                     real64 (& dG)[8] ) const override final
  {

    // The plastic potential function is: G = invQ + invP * dilation

    real64 dG_dP = m_dilation[k];
    real64 dG_dQ = 1.0;

    real64 dG_dP_dP = 0.0;
    real64 dG_dP_dQ = 0.0;
    real64 dG_dP_dH = 0.0;

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
                            localIndex const q,
                            real64 const state ) const override final
  {

    // The hardening function is: cohesion = m_cohesion[k][q] + state * hardeningRate

    return m_cohesion[k][q] + state * m_hardening[k];
  }

  virtual real64 hardeningDerivatives( localIndex const k,
                                       localIndex const GEOSX_UNUSED_PARAM( q ),
                                       real64 const GEOSX_UNUSED_PARAM( state ) ) const override final
  {

    // The hardening function is: cohesion = m_cohesion[k][q] + state * hardeningRate

    return m_hardening[k]; 
  }

  virtual real64 tmpState( localIndex const GEOSX_UNUSED_PARAM( k ),
                           localIndex const GEOSX_UNUSED_PARAM( q ),
                           real64 const GEOSX_UNUSED_PARAM( invP ),
                           real64 const GEOSX_UNUSED_PARAM( invQ ),
                           real64 const plasticMultiplier ) const override final
  {
    // The temporal state variable is updated inside the Newton loops by: state += plasticMultiplier
    // starting from the previous saved state

    return plasticMultiplier; 
  }

  virtual real64 stateDerivatives( localIndex const GEOSX_UNUSED_PARAM( k ),
                                   localIndex const GEOSX_UNUSED_PARAM( q ),
                                   real64 const GEOSX_UNUSED_PARAM( invP ),
                                   real64 const GEOSX_UNUSED_PARAM( invQ ),
                                   real64 const GEOSX_UNUSED_PARAM( plasticMultiplier ) ) const override final
  {
    return 1.0;
  }

  virtual real64 getStateVariable( localIndex const k,
                                   localIndex const q ) const override final
  {
    return m_state[k][q];
  }

  virtual void saveStateVariable( localIndex const k,
                                  localIndex const q,
                                  real64 const state ) const override final
  {
    m_state[k][q] = state;
  }

  /// A reference to the ArrayView holding the friction angle for each element.
  arrayView1d< real64 const > const m_friction;

  /// A reference to the ArrayView holding the dilation angle for each element.
  arrayView1d< real64 const > const m_dilation;

  /// A reference to the ArrayView holding the hardening rate for each element.
  arrayView1d< real64 const > const m_hardening;

  /// A reference to the ArrayView holding the cohesion for each integration point
  arrayView2d< real64 > const m_cohesion;

  /// A reference to the ArrayView holding the state variable for each integration point
  arrayView2d< real64 > const m_state;

};

/**
 * @class DruckerPrager
 *
 * Drucker-Prager material model.
 */
class DruckerPrager : public ElasticIsotropic
{
public:

  /// @typedef Alias for DruckerPragerUpdates
  using KernelWrapper = DruckerPragerUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DruckerPrager( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DruckerPrager() override;


  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "DruckerPrager";

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
    /// string/key for default friction angle
    static constexpr auto defaultFrictionAngleString = "defaultFrictionAngle";

    /// string/key for default dilation angle
    static constexpr auto defaultDilationAngleString = "defaultDilationAngle";

    /// string/key for default hardening rate
    static constexpr auto defaultHardeningString = "defaultHardeningRate";

    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";

    /// string/key for friction angle
    static constexpr auto frictionString  = "friction";

    /// string/key for dilation angle
    static constexpr auto dilationString  = "dilation";

    /// string/key for cohesion
    static constexpr auto hardeningString  = "hardening";

    /// string/key for cohesion
    static constexpr auto cohesionString  = "cohesion";

    /// string/key for state variable
    static constexpr auto stateString  = "stateVariable";
  };

  /**
   * @brief Create a instantiation of the DruckerPragerUpdate class that refers to the data in this.
   * @return An instantiation of DruckerPragerUpdate.
   */
  DruckerPragerUpdates createKernelUpdates() const
  {
    return DruckerPragerUpdates( m_friction,
                                 m_dilation,
                                 m_hardening,
                                 m_cohesion,
                                 m_state,
                                 m_bulkModulus,
                                 m_shearModulus,
                                 m_newStress,
                                 m_oldStress );
  }

protected:
  virtual void postProcessInput() override;

  /// Material parameter: The default value of yield surface slope
  real64 m_defaultFrictionAngle;

  /// Material parameter: The default value of plastic potential slope
  real64 m_defaultDilationAngle;

  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;

  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardening;

  /// Material parameter: The yield surface slope for each element
  array1d< real64 > m_friction;

  /// Material parameter: The plastic potential slope for each element
  array1d< real64 > m_dilation;

  /// Material parameter: The hardening rate each element
  array1d< real64 > m_hardening;

  /// State variable: The current cohesion parameter for each quadrature point
  array2d< real64 > m_cohesion;

  /// State variable: The current equivalent plastic multiplier for each quadrature point
  array2d< real64 > m_state;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP_ */
