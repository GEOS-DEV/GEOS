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

  arrayView2d< real64 const > const getOldPorosity() const
  {
    return getBasePorosityModel().getOldPorosity();
  }

  arrayView2d< real64 const > const getPorosity() const
  {
    return getBasePorosityModel().getPorosity();
  }

  arrayView2d< real64 const > const  getDporosity_dPressure() const
  { return getBasePorosityModel().dPorosity_dPressure(); }

protected:

  // the name of the porosity model
  string m_solidModelName;

  // the name of the porosity model
  string m_porosityModelName;

  // the name of the porosity model
  string m_permeabilityModelName;

private:

  PorosityBase const & getBasePorosityModel() const
  { return this->getParent().template getGroup< PorosityBase >( m_porosityModelName ); }

  PermeabilityBase const & getBasePermModel() const
  { return this->getParent().template getGroup< PermeabilityBase >( m_permeabilityModelName ); }

};

}

}

#endif /* GEOSX_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_ */
