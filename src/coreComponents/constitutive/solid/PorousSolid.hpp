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
 * @file PorousSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE type of the porosity model
 */
template< typename SOLID_TYPE >
class PorousSolidUpdates : public CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  using DiscretizationOps = typename SOLID_TYPE::KernelWrapper::DiscretizationOps;

  /**
   * @brief Constructor
   */
  PorousSolidUpdates( SOLID_TYPE const & solidModel,
                      BiotPorosity const & porosityModel,
                      ConstantPermeability const & permModel ):
    CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >( solidModel, porosityModel, permModel )
  {}

  GEOSX_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & initialPressure,
                          real64 const & pressure,
                          real64 const & deltaPressure,
                          real64 const ( &strainIncrement )[6],
                          real64 ( & totalStress )[6],
                          real64 & dPorosity_dPressure,
                          real64 & dPorosity_dVolStrain,
                          real64 ( & dTotalStress_dPressure )[6],
                          DiscretizationOps & stiffness ) const
  {
    // Compute effective stress and store in totalStress
    m_solidUpdate.smallStrainUpdate( k, q, strainIncrement, totalStress, stiffness );

    updateBiotCoefficient( k );

    // Compute  total stress
    real64 const biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );
    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -biotCoefficient * ( pressure + deltaPressure - initialPressure ) );

    dTotalStress_dPressure[0] = -biotCoefficient;
    dTotalStress_dPressure[1] = -biotCoefficient;
    dTotalStress_dPressure[2] = -biotCoefficient;
    dTotalStress_dPressure[3] = 0;
    dTotalStress_dPressure[4] = 0;
    dTotalStress_dPressure[5] = 0;

    m_porosityUpdate.updateFromPressureAndStrain( k,
                                                  q,
                                                  deltaPressure,
                                                  strainIncrement,
                                                  dPorosity_dPressure,
                                                  dPorosity_dVolStrain );

// TODO uncomment once we start using permeability model in flow.
//    m_permUpdate.updateFromPressureStrain( k,
//                                           q,
//                                           pressure,
//                                           volStrain );
  }

  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k ) const
  {
    // This call is not general like this.
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );

    m_porosityUpdate.updateBiotCoefficient( k, bulkModulus );
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @note If the material model has a strain-dependent material stiffness (e.g.
   * any plasticity, damage, or nonlinear elastic model) then this interface will
   * not work.  Users should instead use one of the interfaces where a strain
   * tensor is provided as input.
   *
   * @param k the element number
   * @param stiffness the stiffness array
   */
  GEOSX_HOST_DEVICE
  void getElasticStiffness( localIndex const k, localIndex const q, real64 ( & stiffness )[6][6] ) const
  {
    m_solidUpdate.getElasticStiffness( k, q, stiffness );
  }

private:

  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_solidUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_porosityUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_permUpdate;

};

/**
 * @brief PorousSolidBase class used for dispatch of all Porous solids.
 */
class PorousSolidBase
{};

/**
 * @brief Class to represent a porous material for poromechanics simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material.
 *
 * @tparam SOLID_TYPE type of solid model
 */
template< typename SOLID_TYPE >
class PorousSolid : public CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = PorousSolidUpdates< SOLID_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PorousSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Porous" ) + SOLID_TYPE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the PorousSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousSolidUpdates.
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
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getSolidModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPorosityModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPermModel;
};



}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_ */
