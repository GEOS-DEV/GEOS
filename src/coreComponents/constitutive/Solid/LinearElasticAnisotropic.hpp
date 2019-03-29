/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 *  @file LinearElasticAnisotropic.hpp
 */

#ifndef LINEARELASTICANISOTROPIC_HPP_
#define LINEARELASTICANISOTROPIC_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const linearElasticIsotropic = "LinearElasticAnisotropic";
}
}

namespace constitutive
{

/**
 * @class LinearElasticAnisotropic
 *
 * Class to provide a linear elastic isotropic material response.
 */
class LinearElasticAnisotropic : public SolidBase
{
public:
  LinearElasticAnisotropic( std::string const & name, ManagedGroup * const parent );

  virtual ~LinearElasticAnisotropic() override;

  virtual std::unique_ptr<ConstitutiveBase>
  DeliverClone( string const & name, ManagedGroup * const parent ) const override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static std::string CatalogName() { return dataRepository::keys::linearElasticIsotropic; }
  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdatePoint( localIndex const k,
                                 localIndex const q,
                                 R2SymTensor const & D,
                                 R2Tensor const & Rot,
                                 integer const systemAssembleFlag ) override;


  void GetStiffness( localIndex const k, real64 c[6][6] ) const;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto stiffness0String  = "stiffness0";
    static constexpr auto stiffnessString  = "stiffness";

  } m_linearElasticIsotropicViewKeys;


protected:
  virtual void PostProcessInput() override;

private:


  SolidBase::StiffnessTensor m_stiffness0;
  array2d<real64> m_stiffness;
};


}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
