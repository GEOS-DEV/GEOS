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
 * @file PorosityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

class PorosityBase : public ConstitutiveBase
{
public:
  PorosityBase( string const & name, Group * const parent );

  virtual ~PorosityBase() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PorosityBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * porosityString() { return "porosity"; }
    static constexpr char const * porosityOldString() { return "porosityOld"; }
    static constexpr char const * dPorosity_dPressureString() { return "dPorosity_dPressure"; }

  } viewKeys;

  arrayView2d< real64 const > const & getPorosity() { return m_porosity.toViewConst(); }

  arrayView2d< real64 const > const & getPorosityOld() { return m_porosityOld.toViewConst(); }

  arrayView2d< real64 > const & getPorosityOld() { return m_porosityOld.toView(); }

  arrayView2d< real64 const > const & dPorosity_dPressure() { return m_dPorosity_dPressure.toViewConst(); }

protected:
  virtual void postProcessInput() override;

  array2d< real64 > m_porosity;

  array2d< real64 > m_porosityOld;

  array2d< real64 > m_dPorosity_dPressure;


};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
