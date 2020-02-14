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
 *  @file LinearElasticAnisotropic.cpp
 */

#include "LinearElasticAnisotropic.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{



LinearElasticAnisotropic::LinearElasticAnisotropic( std::string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultStiffness{}
{
  // These are temporary until we figure out how to read in multidim arrays from input.
  registerWrapper( viewKeyStruct::defaultC11, &(m_defaultStiffness.m_data[0][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 11 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC12, &(m_defaultStiffness.m_data[0][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 12 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC13, &(m_defaultStiffness.m_data[0][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 13 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC14, &(m_defaultStiffness.m_data[0][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 14 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC15, &(m_defaultStiffness.m_data[0][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 15 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC16, &(m_defaultStiffness.m_data[0][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 16 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::defaultC21, &(m_defaultStiffness.m_data[1][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 21 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC22, &(m_defaultStiffness.m_data[1][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 22 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC23, &(m_defaultStiffness.m_data[1][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 23 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC24, &(m_defaultStiffness.m_data[1][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 24 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC25, &(m_defaultStiffness.m_data[1][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 25 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC26, &(m_defaultStiffness.m_data[1][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 26 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::defaultC31, &(m_defaultStiffness.m_data[2][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 31 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC32, &(m_defaultStiffness.m_data[2][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 32 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC33, &(m_defaultStiffness.m_data[2][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 33 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC34, &(m_defaultStiffness.m_data[2][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 34 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC35, &(m_defaultStiffness.m_data[2][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 35 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC36, &(m_defaultStiffness.m_data[2][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 36 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::defaultC41, &(m_defaultStiffness.m_data[3][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 41 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC42, &(m_defaultStiffness.m_data[3][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 42 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC43, &(m_defaultStiffness.m_data[3][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 43 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC44, &(m_defaultStiffness.m_data[3][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 44 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC45, &(m_defaultStiffness.m_data[3][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 45 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC46, &(m_defaultStiffness.m_data[3][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 46 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::defaultC51, &(m_defaultStiffness.m_data[4][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 51 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC52, &(m_defaultStiffness.m_data[4][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 52 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC53, &(m_defaultStiffness.m_data[4][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 53 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC54, &(m_defaultStiffness.m_data[4][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 54 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC55, &(m_defaultStiffness.m_data[4][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 55 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC56, &(m_defaultStiffness.m_data[4][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 56 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::defaultC61, &(m_defaultStiffness.m_data[5][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 61 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC62, &(m_defaultStiffness.m_data[5][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 62 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC63, &(m_defaultStiffness.m_data[5][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 63 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC64, &(m_defaultStiffness.m_data[5][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 64 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC65, &(m_defaultStiffness.m_data[5][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 65 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::defaultC66, &(m_defaultStiffness.m_data[5][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 66 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::defaultStiffnessString, &m_defaultStiffness, 0 )->
//    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default Elastic Stiffness Tensor in Voigt notation");


  // These are temporary until we figure out how to read in multidim arrays from input.
  registerWrapper( viewKeyStruct::c11, &(m_c00), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 11 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c12, &(m_c01), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 12 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c13, &(m_c02), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 13 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c14, &(m_c03), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 14 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c15, &(m_c04), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 15 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c16, &(m_c05), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 16 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::c21, &(m_c10), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 21 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c22, &(m_c11), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 22 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c23, &(m_c12), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 23 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c24, &(m_c13), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 24 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c25, &(m_c14), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 25 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c26, &(m_c15), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 26 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::c31, &(m_c20), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 31 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c32, &(m_c21), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 32 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c33, &(m_c22), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 33 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c34, &(m_c23), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 34 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c35, &(m_c24), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 35 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c36, &(m_c25), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 36 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::c41, &(m_c30), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 41 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c42, &(m_c31), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 42 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c43, &(m_c32), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 43 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c44, &(m_c33), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 44 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c45, &(m_c34), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 45 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c46, &(m_c35), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 46 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::c51, &(m_c40), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 51 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c52, &(m_c41), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 52 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c53, &(m_c42), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 53 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c54, &(m_c43), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 54 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c55, &(m_c44), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 55 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c56, &(m_c45), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 56 component of the Elastic Stiffness Tensor in Voigt notation");

  registerWrapper( viewKeyStruct::c61, &(m_c50), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 61 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c62, &(m_c51), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 62 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c63, &(m_c52), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 63 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c64, &(m_c53), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 64 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c65, &(m_c54), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 65 component of the Elastic Stiffness Tensor in Voigt notation");
  registerWrapper( viewKeyStruct::c66, &(m_c55), 0 )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("The 66 component of the Elastic Stiffness Tensor in Voigt notation");

//  registerWrapper( viewKeyStruct::stiffnessString, &m_stiffness, 0 )->
////    setApplyDefaultValue(-1)->
//    setDescription("Elastic Stiffness Field in Voigt notation");

}


LinearElasticAnisotropic::~LinearElasticAnisotropic()
{}


void
LinearElasticAnisotropic::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<LinearElasticAnisotropic>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  LinearElasticAnisotropic * const newConstitutiveRelation = dynamic_cast<LinearElasticAnisotropic *>(clone.get());


  newConstitutiveRelation->m_defaultStiffness = m_defaultStiffness;
  //newConstitutiveRelation->m_stiffness = m_stiffness;
}

void LinearElasticAnisotropic::AllocateConstitutiveData( dataRepository::Group * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  //m_stiffness.resize( parent->size() );
  m_c00 = m_defaultStiffness.m_data[0][0];
  m_c01 = m_defaultStiffness.m_data[0][1];
  m_c02 = m_defaultStiffness.m_data[0][2];
  m_c03 = m_defaultStiffness.m_data[0][3];
  m_c04 = m_defaultStiffness.m_data[0][4];
  m_c05 = m_defaultStiffness.m_data[0][5];

  m_c10 = m_defaultStiffness.m_data[1][0];
  m_c11 = m_defaultStiffness.m_data[1][1];
  m_c12 = m_defaultStiffness.m_data[1][2];
  m_c13 = m_defaultStiffness.m_data[1][3];
  m_c14 = m_defaultStiffness.m_data[1][4];
  m_c15 = m_defaultStiffness.m_data[1][5];

  m_c20 = m_defaultStiffness.m_data[2][0];
  m_c21 = m_defaultStiffness.m_data[2][1];
  m_c22 = m_defaultStiffness.m_data[2][2];
  m_c23 = m_defaultStiffness.m_data[2][3];
  m_c24 = m_defaultStiffness.m_data[2][4];
  m_c25 = m_defaultStiffness.m_data[2][5];

  m_c30 = m_defaultStiffness.m_data[3][0];
  m_c31 = m_defaultStiffness.m_data[3][1];
  m_c32 = m_defaultStiffness.m_data[3][2];
  m_c33 = m_defaultStiffness.m_data[3][3];
  m_c34 = m_defaultStiffness.m_data[3][4];
  m_c35 = m_defaultStiffness.m_data[3][5];

  m_c40 = m_defaultStiffness.m_data[4][0];
  m_c41 = m_defaultStiffness.m_data[4][1];
  m_c42 = m_defaultStiffness.m_data[4][2];
  m_c43 = m_defaultStiffness.m_data[4][3];
  m_c44 = m_defaultStiffness.m_data[4][4];
  m_c45 = m_defaultStiffness.m_data[4][5];

  m_c50 = m_defaultStiffness.m_data[5][0];
  m_c51 = m_defaultStiffness.m_data[5][1];
  m_c52 = m_defaultStiffness.m_data[5][2];
  m_c53 = m_defaultStiffness.m_data[5][3];
  m_c54 = m_defaultStiffness.m_data[5][4];
  m_c55 = m_defaultStiffness.m_data[5][5];
}

void LinearElasticAnisotropic::PostProcessInput()
{
  getWrapper<array1d<real64>>(viewKeyStruct::c11)->setDefaultValue(m_defaultStiffness.m_data[0][0]);
  getWrapper<array1d<real64>>(viewKeyStruct::c12)->setDefaultValue(m_defaultStiffness.m_data[0][1]);
  getWrapper<array1d<real64>>(viewKeyStruct::c13)->setDefaultValue(m_defaultStiffness.m_data[0][2]);
  getWrapper<array1d<real64>>(viewKeyStruct::c14)->setDefaultValue(m_defaultStiffness.m_data[0][3]);
  getWrapper<array1d<real64>>(viewKeyStruct::c15)->setDefaultValue(m_defaultStiffness.m_data[0][4]);
  getWrapper<array1d<real64>>(viewKeyStruct::c16)->setDefaultValue(m_defaultStiffness.m_data[0][5]);

  getWrapper<array1d<real64>>(viewKeyStruct::c21)->setDefaultValue(m_defaultStiffness.m_data[1][0]);
  getWrapper<array1d<real64>>(viewKeyStruct::c22)->setDefaultValue(m_defaultStiffness.m_data[1][1]);
  getWrapper<array1d<real64>>(viewKeyStruct::c23)->setDefaultValue(m_defaultStiffness.m_data[1][2]);
  getWrapper<array1d<real64>>(viewKeyStruct::c24)->setDefaultValue(m_defaultStiffness.m_data[1][3]);
  getWrapper<array1d<real64>>(viewKeyStruct::c25)->setDefaultValue(m_defaultStiffness.m_data[1][4]);
  getWrapper<array1d<real64>>(viewKeyStruct::c26)->setDefaultValue(m_defaultStiffness.m_data[1][5]);

  getWrapper<array1d<real64>>(viewKeyStruct::c31)->setDefaultValue(m_defaultStiffness.m_data[2][0]);
  getWrapper<array1d<real64>>(viewKeyStruct::c32)->setDefaultValue(m_defaultStiffness.m_data[2][1]);
  getWrapper<array1d<real64>>(viewKeyStruct::c33)->setDefaultValue(m_defaultStiffness.m_data[2][2]);
  getWrapper<array1d<real64>>(viewKeyStruct::c34)->setDefaultValue(m_defaultStiffness.m_data[2][3]);
  getWrapper<array1d<real64>>(viewKeyStruct::c35)->setDefaultValue(m_defaultStiffness.m_data[2][4]);
  getWrapper<array1d<real64>>(viewKeyStruct::c36)->setDefaultValue(m_defaultStiffness.m_data[2][5]);

  getWrapper<array1d<real64>>(viewKeyStruct::c41)->setDefaultValue(m_defaultStiffness.m_data[3][0]);
  getWrapper<array1d<real64>>(viewKeyStruct::c42)->setDefaultValue(m_defaultStiffness.m_data[3][1]);
  getWrapper<array1d<real64>>(viewKeyStruct::c43)->setDefaultValue(m_defaultStiffness.m_data[3][2]);
  getWrapper<array1d<real64>>(viewKeyStruct::c44)->setDefaultValue(m_defaultStiffness.m_data[3][3]);
  getWrapper<array1d<real64>>(viewKeyStruct::c45)->setDefaultValue(m_defaultStiffness.m_data[3][4]);
  getWrapper<array1d<real64>>(viewKeyStruct::c46)->setDefaultValue(m_defaultStiffness.m_data[3][5]);

  getWrapper<array1d<real64>>(viewKeyStruct::c51)->setDefaultValue(m_defaultStiffness.m_data[4][0]);
  getWrapper<array1d<real64>>(viewKeyStruct::c52)->setDefaultValue(m_defaultStiffness.m_data[4][1]);
  getWrapper<array1d<real64>>(viewKeyStruct::c53)->setDefaultValue(m_defaultStiffness.m_data[4][2]);
  getWrapper<array1d<real64>>(viewKeyStruct::c54)->setDefaultValue(m_defaultStiffness.m_data[4][3]);
  getWrapper<array1d<real64>>(viewKeyStruct::c55)->setDefaultValue(m_defaultStiffness.m_data[4][4]);
  getWrapper<array1d<real64>>(viewKeyStruct::c56)->setDefaultValue(m_defaultStiffness.m_data[4][5]);

  getWrapper<array1d<real64>>(viewKeyStruct::c61)->setDefaultValue(m_defaultStiffness.m_data[5][0]);
  getWrapper<array1d<real64>>(viewKeyStruct::c62)->setDefaultValue(m_defaultStiffness.m_data[5][1]);
  getWrapper<array1d<real64>>(viewKeyStruct::c63)->setDefaultValue(m_defaultStiffness.m_data[5][2]);
  getWrapper<array1d<real64>>(viewKeyStruct::c64)->setDefaultValue(m_defaultStiffness.m_data[5][3]);
  getWrapper<array1d<real64>>(viewKeyStruct::c65)->setDefaultValue(m_defaultStiffness.m_data[5][4]);
  getWrapper<array1d<real64>>(viewKeyStruct::c66)->setDefaultValue(m_defaultStiffness.m_data[5][5]);

}

void LinearElasticAnisotropic::StateUpdatePoint( localIndex const k,
                                                 localIndex const q,
                                                 R2SymTensor const & D,
                                                 R2Tensor const & Rot,
                                                 integer const GEOSX_UNUSED_PARAM( updateStiffnessFlag ) )
{
  R2SymTensor T;
  real64 * const GEOSX_RESTRICT Tdata = T.Data();
  real64 const * const GEOSX_RESTRICT Ddata = D.Data();
//  real64 const (&c)[6][6] = m_stiffness[k].m_data;

  Tdata[0] += m_c00[k]*Ddata[0] + m_c01[k]*Ddata[2] + m_c02[k]*Ddata[5] + m_c03[k]*(2*Ddata[4]) + m_c04[k]*(2*Ddata[3]) + m_c05[k]*(2*Ddata[1]);
  Tdata[2] += m_c10[k]*Ddata[0] + m_c11[k]*Ddata[2] + m_c12[k]*Ddata[5] + m_c13[k]*(2*Ddata[4]) + m_c14[k]*(2*Ddata[3]) + m_c15[k]*(2*Ddata[1]);
  Tdata[5] += m_c20[k]*Ddata[0] + m_c21[k]*Ddata[2] + m_c22[k]*Ddata[5] + m_c23[k]*(2*Ddata[4]) + m_c24[k]*(2*Ddata[3]) + m_c25[k]*(2*Ddata[1]);
  Tdata[4] += m_c30[k]*Ddata[0] + m_c31[k]*Ddata[2] + m_c32[k]*Ddata[5] + m_c33[k]*(2*Ddata[4]) + m_c34[k]*(2*Ddata[3]) + m_c35[k]*(2*Ddata[1]);
  Tdata[3] += m_c40[k]*Ddata[0] + m_c41[k]*Ddata[2] + m_c42[k]*Ddata[5] + m_c43[k]*(2*Ddata[4]) + m_c44[k]*(2*Ddata[3]) + m_c45[k]*(2*Ddata[1]);
  Tdata[1] += m_c50[k]*Ddata[0] + m_c51[k]*Ddata[2] + m_c52[k]*Ddata[5] + m_c53[k]*(2*Ddata[4]) + m_c54[k]*(2*Ddata[3]) + m_c55[k]*(2*Ddata[1]);

  R2SymTensor temp2 = m_stress[ k ][ q ];
  temp2 += T;

  T.QijAjkQlk( temp2, Rot );
  
  for ( localIndex i = 0; i < 6; ++i )
  {
    m_stress( k, q, i ) = T.Data()[ i ];
  }
}

//void LinearElasticAnisotropic::GetStiffness( localIndex const k, real64 c[6][6] ) const
//{
//  for( int i=0 ; i<6 ; ++i )
//  {
//    for( int j=0 ; j<6 ; ++j )
//    {
//      c[i][j] = m_stiffness[k].m_data[i][j];
//    }
//  }
//}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticAnisotropic, std::string const &, Group * const )
}
} /* namespace geosx */
