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
 * @file PermeabilityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

class PermeabilityBase : public ConstitutiveBase
{
public:
  PermeabilityBase( string const & name, Group * const parent );

  virtual ~PermeabilityBase() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PermeabilityBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto permeabilityString = "permeability";
  } viewKeys;

protected:
  virtual void postProcessInput() override;

  array3d< real64 > m_permeability;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_
