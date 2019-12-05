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

#ifndef GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICANISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICANISOTROPIC_HPP_
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
  LinearElasticAnisotropic( string const & name, Group * const parent );

  virtual ~LinearElasticAnisotropic() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr<ConstitutiveBase> & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr auto m_catalogNameString = "LinearElasticAnisotropic";
  static std::string CatalogName() { return m_catalogNameString; }
  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdatePoint( localIndex const k,
                                 localIndex const q,
                                 R2SymTensor const & D,
                                 R2Tensor const & Rot,
                                 integer const updateStiffnessFlag ) override;


//  void GetStiffness( localIndex const k, real64 c[6][6] ) const;

  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    static constexpr auto defaultC11 = "defaultC11";
    static constexpr auto defaultC12 = "defaultC12";
    static constexpr auto defaultC13 = "defaultC13";
    static constexpr auto defaultC14 = "defaultC14";
    static constexpr auto defaultC15 = "defaultC15";
    static constexpr auto defaultC16 = "defaultC16";

    static constexpr auto defaultC21 = "defaultC21";
    static constexpr auto defaultC22 = "defaultC22";
    static constexpr auto defaultC23 = "defaultC23";
    static constexpr auto defaultC24 = "defaultC24";
    static constexpr auto defaultC25 = "defaultC25";
    static constexpr auto defaultC26 = "defaultC26";

    static constexpr auto defaultC31 = "defaultC31";
    static constexpr auto defaultC32 = "defaultC32";
    static constexpr auto defaultC33 = "defaultC33";
    static constexpr auto defaultC34 = "defaultC34";
    static constexpr auto defaultC35 = "defaultC35";
    static constexpr auto defaultC36 = "defaultC36";

    static constexpr auto defaultC41 = "defaultC41";
    static constexpr auto defaultC42 = "defaultC42";
    static constexpr auto defaultC43 = "defaultC43";
    static constexpr auto defaultC44 = "defaultC44";
    static constexpr auto defaultC45 = "defaultC45";
    static constexpr auto defaultC46 = "defaultC46";

    static constexpr auto defaultC51 = "defaultC51";
    static constexpr auto defaultC52 = "defaultC52";
    static constexpr auto defaultC53 = "defaultC53";
    static constexpr auto defaultC54 = "defaultC54";
    static constexpr auto defaultC55 = "defaultC55";
    static constexpr auto defaultC56 = "defaultC56";

    static constexpr auto defaultC61 = "defaultC61";
    static constexpr auto defaultC62 = "defaultC62";
    static constexpr auto defaultC63 = "defaultC63";
    static constexpr auto defaultC64 = "defaultC64";
    static constexpr auto defaultC65 = "defaultC65";
    static constexpr auto defaultC66 = "defaultC66";


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

  /**
   * @struct wrapper struct for the stiffness tensor
   */
  struct StiffnessTensor
  {
    real64 m_data[6][6];
  };

//  /**
//   * @brief accessor for the stiffness field
//   * @return
//   */
//  arrayView1d<StiffnessTensor> const &       stiffness()       { return m_stiffness; }
//
//  /**
//   * @brief immutable accessor for the stiffness field
//   * @return
//   */
//  arrayView1d<StiffnessTensor const> const & stiffness() const { return m_stiffness; }


  void setDefaultStiffness( StiffnessTensor const & input )
  {
    m_defaultStiffness = input;
  }

  class KernelWrapper
  {
  public:
    KernelWrapper( arrayView1d<real64 const> const & c00,
                   arrayView1d<real64 const> const & c01,
                   arrayView1d<real64 const> const & c02,
                   arrayView1d<real64 const> const & c03,
                   arrayView1d<real64 const> const & c04,
                   arrayView1d<real64 const> const & c05,
                   arrayView1d<real64 const> const & c10,
                   arrayView1d<real64 const> const & c11,
                   arrayView1d<real64 const> const & c12,
                   arrayView1d<real64 const> const & c13,
                   arrayView1d<real64 const> const & c14,
                   arrayView1d<real64 const> const & c15,
                   arrayView1d<real64 const> const & c20,
                   arrayView1d<real64 const> const & c21,
                   arrayView1d<real64 const> const & c22,
                   arrayView1d<real64 const> const & c23,
                   arrayView1d<real64 const> const & c24,
                   arrayView1d<real64 const> const & c25,
                   arrayView1d<real64 const> const & c30,
                   arrayView1d<real64 const> const & c31,
                   arrayView1d<real64 const> const & c32,
                   arrayView1d<real64 const> const & c33,
                   arrayView1d<real64 const> const & c34,
                   arrayView1d<real64 const> const & c35,
                   arrayView1d<real64 const> const & c40,
                   arrayView1d<real64 const> const & c41,
                   arrayView1d<real64 const> const & c42,
                   arrayView1d<real64 const> const & c43,
                   arrayView1d<real64 const> const & c44,
                   arrayView1d<real64 const> const & c45,
                   arrayView1d<real64 const> const & c50,
                   arrayView1d<real64 const> const & c51,
                   arrayView1d<real64 const> const & c52,
                   arrayView1d<real64 const> const & c53,
                   arrayView1d<real64 const> const & c54,
                   arrayView1d<real64 const> const & c55 ) :
      m_stiffnessView{ {c00,c01,c02,c03,c04,c05},
                       {c10,c11,c12,c13,c14,c15},
                       {c20,c21,c22,c23,c24,c25},
                       {c30,c31,c32,c33,c34,c35},
                       {c40,c41,c42,c43,c44,c45},
                       {c50,c51,c52,c53,c54,c55} }
    {}

    /**
     * accessor to return the stiffness at a given element
     * @param k the element number
     * @param c the stiffness array
     */
    GEOSX_HOST_DEVICE inline
    void GetStiffness( localIndex const k, real64 (&c)[6][6] ) const
    {
      for( int i=0 ; i<6 ; ++i )
      {
        for( int j=0 ; j<6 ; ++j )
        {
          c[i][j] = m_stiffnessView[i][j](k);
        }
      }
    }

  private:
    arrayView1d<real64 const> const m_stiffnessView[6][6];

  };

  KernelWrapper createKernelWrapper() const
  { return KernelWrapper( m_c00, m_c01, m_c02, m_c03, m_c04, m_c05,
                          m_c10, m_c11, m_c12, m_c13, m_c14, m_c15,
                          m_c20, m_c21, m_c22, m_c23, m_c24, m_c25,
                          m_c30, m_c31, m_c32, m_c33, m_c34, m_c35,
                          m_c40, m_c41, m_c42, m_c43, m_c44, m_c45,
                          m_c50, m_c51, m_c52, m_c53, m_c54, m_c55 ); }
protected:
  virtual void PostProcessInput() override;

private:


  StiffnessTensor m_defaultStiffness; /// default value for stiffness tensor
//  array1d<StiffnessTensor> m_stiffness; /// stiffness tensor field

  array1d<real64> m_c00;
  array1d<real64> m_c01;
  array1d<real64> m_c02;
  array1d<real64> m_c03;
  array1d<real64> m_c04;
  array1d<real64> m_c05;

  array1d<real64> m_c10;
  array1d<real64> m_c11;
  array1d<real64> m_c12;
  array1d<real64> m_c13;
  array1d<real64> m_c14;
  array1d<real64> m_c15;

  array1d<real64> m_c20;
  array1d<real64> m_c21;
  array1d<real64> m_c22;
  array1d<real64> m_c23;
  array1d<real64> m_c24;
  array1d<real64> m_c25;

  array1d<real64> m_c30;
  array1d<real64> m_c31;
  array1d<real64> m_c32;
  array1d<real64> m_c33;
  array1d<real64> m_c34;
  array1d<real64> m_c35;

  array1d<real64> m_c40;
  array1d<real64> m_c41;
  array1d<real64> m_c42;
  array1d<real64> m_c43;
  array1d<real64> m_c44;
  array1d<real64> m_c45;

  array1d<real64> m_c50;
  array1d<real64> m_c51;
  array1d<real64> m_c52;
  array1d<real64> m_c53;
  array1d<real64> m_c54;
  array1d<real64> m_c55;


};


}

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_LINEARELASTICANISOTROPIC_HPP_ */
