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
 *  @file LinearElasticIsotropic.hpp
 */

#ifndef LINEARELASTICISOTROPIC_HPP_
#define LINEARELASTICISOTROPIC_HPP_
#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class LinearElasticIsotropic
 *
 * Class to provide a linear elastic isotropic material response.
 */
class LinearElasticIsotropic : public SolidBase
{
public:
  LinearElasticIsotropic( string const & name, ManagedGroup * const parent );

  virtual ~LinearElasticIsotropic() override;

  virtual void
  DeliverClone( string const & name,
                ManagedGroup * const parent,
                std::unique_ptr<ConstitutiveBase> & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr auto m_catalogNameString = "LinearElasticIsotropic";
  static std::string CatalogName() { return m_catalogNameString; }
  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdatePoint( localIndex const k,
                                 localIndex const q,
                                 R2SymTensor const & D,
                                 R2Tensor const & Rot,
                                 integer const systemAssembleFlag ) override;


  void GetStiffness( localIndex const k, real64 c[6][6] ) const;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto bulkModulus0String  = "BulkModulus0";
    static constexpr auto poissonRatioString =  "PoissonRatio" ;
    static constexpr auto shearModulus0String = "ShearModulus0";
    static constexpr auto youngsModulus0String =  "YoungsModulus" ;

    static constexpr auto compressibilityString =  "compressibility" ;
    static constexpr auto referencePressureString =  "referencePressure" ;
    static constexpr auto biotCoefficientString =  "BiotCoefficient" ;


    static constexpr auto bulkModulusString  = "BulkModulus";
    static constexpr auto shearModulusString = "ShearModulus";


  } m_linearElasticIsotropicViewKeys;


protected:
  virtual void PostProcessInput() override;

private:


  real64 m_bulkModulus0;
  real64 m_shearModulus0;
  array1d<real64> m_bulkModulus;
  array1d<real64> m_shearModulus;

};


}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
