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
 * @file PorousElasticTransverseIsotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROUSELASTICTRANSVERSEISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROUSELASTICTRANSVERSEISOTROPIC_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/ElasticTransverseIsotropic.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 */
class PorousElasticTransverseIsotropicUpdates : public PorousSolidUpdates< ElasticTransverseIsotropic >
{
public:

  /**
   * @brief Constructor
   */
  PorousElasticTransverseIsotropicUpdates( ElasticTransverseIsotropic const & solidModel,
                                           BiotPorosity const & porosityModel,
                                           ConstantPermeability const & permModel ):
    PorousSolidUpdates< ElasticTransverseIsotropic >( solidModel, porosityModel, permModel )
  {}

  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k ) const
  {
    real64 const c11 = m_solidUpdate.getC11( k );
    real64 const c13 = m_solidUpdate.getC11( k );
    real64 const c33 = m_solidUpdate.getC11( k );
    real64 const c66 = m_solidUpdate.getC11( k );

    real64 const c12 = c11 - 2.0 * c66;

    m_porosityUpdate.updateBiotCoefficient( k, c11, c12, c13, c33 );
  }
};

/**
 * @brief PorousElasticTransverseIsotropicBase class used for dispatch of all Porous solids.
 */
class PorousElasticTransverseIsotropicBase : public SolidBase
{};

/**
 * @brief Class to represent a porous material for poromechanics simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material.
 *
 * @tparam ElasticTransverseIsotropic type of solid model
 */
class PorousElasticTransverseIsotropic : public CoupledSolid< ElasticTransverseIsotropic, BiotPorosity, ConstantPermeability >
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = PorousElasticTransverseIsotropicUpdates;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousElasticTransverseIsotropic( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PorousElasticTransverseIsotropic() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Porous" ) + ElasticTransverseIsotropic::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the PorousElasticTransverseIsotropicUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousElasticTransverseIsotropicUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( getSolidModel(),
                          getPorosityModel(),
                          getPermModel() );
  }

  /**
   * @brief Const/non-mutable accessor for density
   * @return Accessor
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return getSolidModel().getDensity();
  }

private:

  using CoupledSolid< ElasticTransverseIsotropic, BiotPorosity, ConstantPermeability >::getSolidModel;
  using CoupledSolid< ElasticTransverseIsotropic, BiotPorosity, ConstantPermeability >::getPorosityModel;
  using CoupledSolid< ElasticTransverseIsotropic, BiotPorosity, ConstantPermeability >::getPermModel;
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROUSELASTICTRANSVERSEISOTROPIC_HPP_ */
