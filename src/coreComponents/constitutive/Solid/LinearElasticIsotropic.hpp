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
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const linearElasticIsotropic = "LinearElasticIsotropic";
}
}

namespace constitutive
{

class LinearElasticIsotropic : public ConstitutiveBase
{
public:
  LinearElasticIsotropic( std::string const & name, ManagedGroup * const parent );

  virtual ~LinearElasticIsotropic() override;

  virtual std::unique_ptr<ConstitutiveBase>
  DeliverClone( string const & name, ManagedGroup * const parent ) const override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static std::string CatalogName() { return dataRepository::keys::linearElasticIsotropic; }
  virtual string GetCatalogName() override { return CatalogName(); }

  R2SymTensor  StateUpdatePoint( R2SymTensor const & D,
                                 R2Tensor const & Rot,
                                 localIndex const i,
                                 localIndex const q,
                                 integer const systemAssembleFlag ) override;


  void GetStiffness( localIndex const k, real64 c[6][6] ) const;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto bulkModulus0String  = "BulkModulus0";
    static constexpr auto density0String  = "density0";
    static constexpr auto poissonRatioString =  "PoissonRatio" ;
    static constexpr auto shearModulus0String = "ShearModulus0";
    static constexpr auto youngsModulus0String =  "YoungsModulus" ;

    static constexpr auto compressibilityString =  "compressibility" ;
    static constexpr auto referencePressureString =  "referencePressure" ;
    static constexpr auto biotCoefficientString =  "BiotCoefficient" ;


    static constexpr auto bulkModulusString  = "BulkModulus";
    static constexpr auto densityString  = "density";
    static constexpr auto deviatorStressString = "DeviatorStress";
    static constexpr auto meanStressString = "MeanStress";
    static constexpr auto shearModulusString = "ShearModulus";


  } m_linearElasticIsotropicViewKeys;


  virtual void StateUpdatePointPressure(real64 const & pres,
                                        localIndex const k,
                                        localIndex const q) override final;

protected:
  virtual void PostProcessInput() override;

private:


  real64 m_bulkModulus0;
  real64 m_shearModulus0;
  real64 m_density0;
  array2d<real64> m_density;
  array1d<real64> m_bulkModulus;
  array1d<real64> m_shearModulus;

  array2d<real64> m_meanStress;
  array2d<R2SymTensor> m_deviatorStress;

  /// scalar compressibility parameter
  real64 m_compressibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// scalar Biot's coefficient
  real64 m_biotCoefficient;

  array2d<real64> m_poreVolumeMultiplier;
  array2d<real64> m_dPVMult_dPressure;

  ExponentialRelation<real64, ExponentApproximationType::Linear> m_poreVolumeRelation;
};

inline void LinearElasticIsotropic::StateUpdatePointPressure( real64 const & pres,
                                                              localIndex const k,
                                                              localIndex const q )
{
  m_poreVolumeRelation.Compute( pres, m_poreVolumeMultiplier[k][q], m_dPVMult_dPressure[k][q] );
}

}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
