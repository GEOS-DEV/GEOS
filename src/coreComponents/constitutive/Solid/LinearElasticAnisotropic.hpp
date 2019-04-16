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


namespace constitutive
{

/**
 * @class LinearElasticAnisotropic
 *
 * Class to provide a linear elastic anisotropic material response.
 */
class LinearElasticAnisotropic : public SolidBase
{
public:
  LinearElasticAnisotropic( string const & name, ManagedGroup * const parent );

  virtual ~LinearElasticAnisotropic() override;

  virtual void
  DeliverClone( string const & name,
                ManagedGroup * const parent,
                std::unique_ptr<ConstitutiveBase> & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr auto m_catalogNameString = "LinearElasticAnisotropic";
  static std::string CatalogName() { return m_catalogNameString; }
  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdatePoint( localIndex const k,
                                 localIndex const q,
                                 R2SymTensor const & D,
                                 R2Tensor const & Rot,
                                 integer const systemAssembleFlag ) override;


  void GetStiffness( localIndex const k, real64 c[6][6] ) const;

  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    static constexpr auto c11 = "c11";
    static constexpr auto c12 = "c12";
    static constexpr auto c13 = "c13";
    static constexpr auto c14 = "c14";
    static constexpr auto c15 = "c15";
    static constexpr auto c16 = "c16";

    static constexpr auto c21 = "c21";
    static constexpr auto c22 = "c22";
    static constexpr auto c23 = "c23";
    static constexpr auto c24 = "c24";
    static constexpr auto c25 = "c25";
    static constexpr auto c26 = "c26";

    static constexpr auto c31 = "c31";
    static constexpr auto c32 = "c32";
    static constexpr auto c33 = "c33";
    static constexpr auto c34 = "c34";
    static constexpr auto c35 = "c35";
    static constexpr auto c36 = "c36";

    static constexpr auto c41 = "c41";
    static constexpr auto c42 = "c42";
    static constexpr auto c43 = "c43";
    static constexpr auto c44 = "c44";
    static constexpr auto c45 = "c45";
    static constexpr auto c46 = "c46";

    static constexpr auto c51 = "c51";
    static constexpr auto c52 = "c52";
    static constexpr auto c53 = "c53";
    static constexpr auto c54 = "c54";
    static constexpr auto c55 = "c55";
    static constexpr auto c56 = "c56";

    static constexpr auto c61 = "c61";
    static constexpr auto c62 = "c62";
    static constexpr auto c63 = "c63";
    static constexpr auto c64 = "c64";
    static constexpr auto c65 = "c65";
    static constexpr auto c66 = "c66";

    static constexpr auto defaultStiffnessString  = "defaultStiffness";
    static constexpr auto stiffnessString  = "stiffness";

  };

  struct StiffnessTensor
  {
    real64 m_data[6][6];
  };

  arrayView1d<StiffnessTensor> const &       stiffness()       { return m_stiffness; }
  arrayView1d<StiffnessTensor const> const & stiffness() const { return m_stiffness; }


protected:
  virtual void PostProcessInput() override;

private:


  StiffnessTensor m_defaultStiffness;
  array1d<StiffnessTensor> m_stiffness;
};


}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
