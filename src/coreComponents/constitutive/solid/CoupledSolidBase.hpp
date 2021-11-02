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
 * @file CoupledSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"

namespace geosx
{
namespace constitutive
{


class CoupledSolidBase : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  CoupledSolidBase( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~CoupledSolidBase() override;

  struct viewKeyStruct
  {
    static constexpr char const * solidModelNameString() { return "solidModelName"; }
    static constexpr char const * porosityModelNameString() { return "porosityModelName"; }
    static constexpr char const * permeabilityModelNameString() { return "permeabilityModelName"; }
  };

  virtual std::vector< string > getSubRelationNames() const override final
  {
    return { m_solidModelName,
             m_porosityModelName,
             m_permeabilityModelName };

  }

  /**
   * @brief get the old porosity.
   * return a constant arrayView2d to the old porosity
   */
  arrayView2d< real64 const > const getOldPorosity() const
  {
    return getBasePorosityModel().getOldPorosity();
  }

  /*
   * @brief get the new porosity.
   * return a constant arrayView2d to the new porosity
   */
  arrayView2d< real64 const > const getPorosity() const
  {
    return getBasePorosityModel().getPorosity();
  }

  /**
   * @brief get the reference porosity.
   * return a constant arrayView1d to the reference porosity
   */
  arrayView1d< real64 const > const getReferencePorosity() const
  {
    return getBasePorosityModel().getReferencePorosity();
  }

  /**
   * @brief get the dPorosity_dPressure.
   * return a constant arrayView2d to dPorosity_dPressure
   */
  arrayView2d< real64 const > const  getDporosity_dPressure() const
  { return getBasePorosityModel().dPorosity_dPressure(); }

  /**
   * @brief initialize the constitutive models fields.
   */
  void initializeState() const
  {
    getBasePorosityModel().initializeState();
  }

  virtual void saveConvergedState() const override final
  {
    getBasePorosityModel().saveConvergedState();
  }
protected:

  /// the name of the solid model
  string m_solidModelName;

  /// the name of the porosity model
  string m_porosityModelName;

  /// the name of the permeability model
  string m_permeabilityModelName;

private:

  /**
   * @brief get a PorosityBase constant reference to the porosity model
   * return a constant PorosityBase reference to the porosity model
   */
  PorosityBase const & getBasePorosityModel() const
  { return this->getParent().template getGroup< PorosityBase >( m_porosityModelName ); }
  /**
   * @brief get a Permeability base constant reference to the permeability model
   * return a constant PermeabilityBase reference to the permeability model
   */
  PermeabilityBase const & getBasePermModel() const
  { return this->getParent().template getGroup< PermeabilityBase >( m_permeabilityModelName ); }

};

}

}

#endif /* GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_ */
