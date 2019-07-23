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
  /**
   * constructor
   * @param name name of the instance in the catalog
   * @param parent the group which contains this instance
   */
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
                                 integer const updateStiffnessFlag ) override;


  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    static constexpr auto bulkModulus0String  = "defaultBulkModulus";
    static constexpr auto poissonRatioString =  "defaultPoissonRatio" ;
    static constexpr auto shearModulus0String = "defaultShearModulus";
    static constexpr auto youngsModulus0String =  "defaultYoungsModulus" ;

    static constexpr auto bulkModulusString  = "BulkModulus";
    static constexpr auto shearModulusString = "ShearModulus";
  };


  real64   bulkModulus0()  const { return m_defaultBulkModulus; }
  real64 & bulkModulus0()        { return m_defaultBulkModulus; }

  real64 defaultShearModulus() const { return m_defaultShearModulus; }
  real64 & defaultShearModulus()     { return m_defaultShearModulus; }

  arrayView1d<real64> const &       bulkModulus()       { return m_bulkModulus; }
  arrayView1d<real64 const> const & bulkModulus() const { return m_bulkModulus; }

  arrayView1d<real64> const &       shearModulus()       { return m_shearModulus; }
  arrayView1d<real64 const> const & shearModulus() const { return m_shearModulus; }

  class KernelWrapper
  {
  public:
    KernelWrapper( arrayView1d<real64 const> const bulkModulus,
                   arrayView1d<real64 const> const shearModulus ) :
      m_bulkModulus( bulkModulus ),
      m_shearModulus( shearModulus )
    {}

    /**
     * accessor to return the stiffness at a given element
     * @param k the element number
     * @param c the stiffness array
     */
    GEOSX_HOST_DEVICE inline
    void GetStiffness( localIndex const k, real64 (&c)[6][6] ) const
    {
      real64 const G = m_shearModulus[k];
      real64 const Lame = m_bulkModulus[k] - 2.0/3.0 * G;

      memset( c, 0, sizeof( c ) );

      c[0][0] = Lame + 2 * G;
      c[0][1] = Lame;
      c[0][2] = Lame;

      c[1][0] = Lame;
      c[1][1] = Lame + 2 * G;
      c[1][2] = Lame;

      c[2][0] = Lame;
      c[2][1] = Lame;
      c[2][2] = Lame + 2 * G;

      c[3][3] = G;

      c[4][4] = G;

      c[5][5] = G;
    }

  private:
    arrayView1d<real64 const> const m_bulkModulus;
    arrayView1d<real64 const> const m_shearModulus;
  };

  KernelWrapper createKernelWrapper() const
  { return KernelWrapper( m_bulkModulus, m_shearModulus ); }

protected:
  virtual void PostProcessInput() override;

private:

  real64 m_defaultBulkModulus;
  real64 m_defaultShearModulus;
  array1d<real64> m_bulkModulus;
  array1d<real64> m_shearModulus;
};


}

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_HYPOELASTICLINEAR_HPP_ */
