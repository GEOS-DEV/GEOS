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
 *  @file LinearElasticAnisotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_LINEARVISCOELASTICANISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_LINEARVISCOELASTICANISOTROPIC_HPP_
#include "LinearElasticAnisotropic.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{


namespace constitutive
{

/**
 * @class LinearViscoElasticAnisotropic
 *
 * Class to provide a linear viscoelastic anisotropic material response.
 */
class LinearViscoElasticAnisotropic : public LinearElasticAnisotropic
{
public:
  LinearViscoElasticAnisotropic( string const & name, Group * const parent );

  virtual ~LinearViscoElasticAnisotropic() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr<ConstitutiveBase> & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr auto m_catalogNameString = "LinearViscoElasticAnisotropic";

  virtual void StateUpdatePoint( localIndex const k,
                                 localIndex const q,
                                 R2SymTensor const & D,
                                 R2Tensor const & Rot,
                                 real64 const dt,
                                 integer const updateStiffnessFlag ) override;

  struct viewKeyStruct : public LinearElasticAnisotropic::viewKeyStruct
  {
    static constexpr auto viscosityString =  "viscosity" ;
  };

private:
  /// scalar viscosity parameter
  real64 m_viscosity;

  array3d<real64, solid::STRESS_PERMUTATION> m_elasticStress;
};


}

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_LINEARVISCOELASTICANISOTROPIC_HPP_ */
