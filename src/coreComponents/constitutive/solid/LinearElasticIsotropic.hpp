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
 *  @file LinearElasticIsotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICISOTROPIC_HPP_
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
  LinearElasticIsotropic( string const & name, Group * const parent );

  virtual ~LinearElasticIsotropic() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr<ConstitutiveBase> & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
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
    static constexpr auto defaultBulkModulusString  = "defaultBulkModulus";
    static constexpr auto defaultPoissonRatioString =  "defaultPoissonRatio" ;
    static constexpr auto defaultShearModulusString = "defaultShearModulus";
    static constexpr auto defaultYoungsModulusString =  "defaultYoungsModulus" ;

    static constexpr auto bulkModulusString  = "BulkModulus";
    static constexpr auto shearModulusString = "ShearModulus";
  };

  void setDefaultBulkModulus (real64 const bulkModulus) {m_defaultBulkModulus = bulkModulus;}
  void setDefaultShearModulus (real64 const shearModulus) {m_defaultShearModulus = shearModulus;}

  arrayView1d<real64> const &       bulkModulus()       { return m_bulkModulus; }
  arrayView1d<real64 const> const & bulkModulus() const { return m_bulkModulus; }

  arrayView1d<real64> const &       shearModulus()       { return m_shearModulus; }
  arrayView1d<real64 const> const & shearModulus() const { return m_shearModulus; }

  class KernelWrapper
  {
  public:
    KernelWrapper( arrayView1d<real64 const> const & bulkModulus,
                   arrayView1d<real64 const> const & shearModulus,
                   arrayView2d<R2SymTensor> const & stress) :
      m_bulkModulus( bulkModulus ),
      m_shearModulus( shearModulus ),
      m_stress(stress)
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

    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void StateUpdatePoint( localIndex const k,
                           localIndex const q,
                           R2SymTensor const & D,
                           R2Tensor const & Rot) const;

    arrayView1d<real64 const> const m_bulkModulus;
    arrayView1d<real64 const> const m_shearModulus;
    arrayView2d<R2SymTensor> const m_stress;

  };

  KernelWrapper createKernelWrapper()
  { return KernelWrapper( m_bulkModulus, m_shearModulus, m_stress ); }


protected:
  virtual void PostProcessInput() override;

private:

  real64 m_defaultBulkModulus;
  real64 m_defaultShearModulus;
  array1d<real64> m_bulkModulus;
  array1d<real64> m_shearModulus;
  bool m_postProcessed = false;
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void LinearElasticIsotropic::KernelWrapper::StateUpdatePoint( localIndex const k,
                                                              localIndex const q,
                                                              R2SymTensor const & D,
                                                              R2Tensor const & Rot ) const
{
  real64 meanStresIncrement = D.Trace();

  R2SymTensor temp = D;
  temp.PlusIdentity( -meanStresIncrement / 3.0 );
  temp *= 2.0 * m_shearModulus[k];
  meanStresIncrement *= m_bulkModulus[k];
  temp.PlusIdentity( meanStresIncrement );

  m_stress[k][q] += temp;

  temp.QijAjkQlk( m_stress[k][q], Rot );
  m_stress[k][q] = temp;
}

}

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICISOTROPIC_HPP_ */
