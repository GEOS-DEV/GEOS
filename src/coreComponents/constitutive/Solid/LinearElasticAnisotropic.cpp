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
 *  @file LinearElasticAnisotropic.cpp
 */

#include "LinearElasticAnisotropic.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{



LinearElasticAnisotropic::LinearElasticAnisotropic( std::string const & name, ManagedGroup * const parent ):
  SolidBase( name, parent ),
  m_defaultStiffness{},
  m_stiffness{}
{
  // These are temporary until we figure out how to read in multidim arrays from input.
  RegisterViewWrapper( viewKeyStruct::c11, &(m_defaultStiffness.m_data[0][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 11 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c12, &(m_defaultStiffness.m_data[0][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 12 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c13, &(m_defaultStiffness.m_data[0][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 13 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c14, &(m_defaultStiffness.m_data[0][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 14 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c15, &(m_defaultStiffness.m_data[0][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 15 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c16, &(m_defaultStiffness.m_data[0][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 16 component of the Elastic Stiffness Tensor in Voigt notation");

  RegisterViewWrapper( viewKeyStruct::c21, &(m_defaultStiffness.m_data[1][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 21 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c22, &(m_defaultStiffness.m_data[1][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 22 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c23, &(m_defaultStiffness.m_data[1][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 23 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c24, &(m_defaultStiffness.m_data[1][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 24 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c25, &(m_defaultStiffness.m_data[1][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 25 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c26, &(m_defaultStiffness.m_data[1][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 26 component of the Elastic Stiffness Tensor in Voigt notation");

  RegisterViewWrapper( viewKeyStruct::c31, &(m_defaultStiffness.m_data[2][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 31 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c32, &(m_defaultStiffness.m_data[2][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 32 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c33, &(m_defaultStiffness.m_data[2][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 33 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c34, &(m_defaultStiffness.m_data[2][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 34 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c35, &(m_defaultStiffness.m_data[2][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 35 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c36, &(m_defaultStiffness.m_data[2][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 36 component of the Elastic Stiffness Tensor in Voigt notation");

  RegisterViewWrapper( viewKeyStruct::c41, &(m_defaultStiffness.m_data[3][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 41 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c42, &(m_defaultStiffness.m_data[3][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 42 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c43, &(m_defaultStiffness.m_data[3][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 43 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c44, &(m_defaultStiffness.m_data[3][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 44 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c45, &(m_defaultStiffness.m_data[3][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 45 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c46, &(m_defaultStiffness.m_data[3][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 46 component of the Elastic Stiffness Tensor in Voigt notation");

  RegisterViewWrapper( viewKeyStruct::c51, &(m_defaultStiffness.m_data[4][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 51 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c52, &(m_defaultStiffness.m_data[4][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 52 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c53, &(m_defaultStiffness.m_data[4][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 53 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c54, &(m_defaultStiffness.m_data[4][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 54 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c55, &(m_defaultStiffness.m_data[4][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 55 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c56, &(m_defaultStiffness.m_data[4][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 56 component of the Elastic Stiffness Tensor in Voigt notation");

  RegisterViewWrapper( viewKeyStruct::c61, &(m_defaultStiffness.m_data[5][0]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 61 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c62, &(m_defaultStiffness.m_data[5][1]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 62 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c63, &(m_defaultStiffness.m_data[5][2]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 63 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c64, &(m_defaultStiffness.m_data[5][3]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 64 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c65, &(m_defaultStiffness.m_data[5][4]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 65 component of the Elastic Stiffness Tensor in Voigt notation");
  RegisterViewWrapper( viewKeyStruct::c66, &(m_defaultStiffness.m_data[5][5]), 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default for the 66 component of the Elastic Stiffness Tensor in Voigt notation");

  RegisterViewWrapper( viewKeyStruct::defaultStiffnessString, &m_defaultStiffness, 0 )->
//    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default Elastic Stiffness Tensor in Voigt notation");

  RegisterViewWrapper( viewKeyStruct::stiffnessString, &m_stiffness, 0 )->
//    setApplyDefaultValue(-1)->
    setDescription("Elastic Stiffness Field in Voigt notation");

}


LinearElasticAnisotropic::~LinearElasticAnisotropic()
{}


void
LinearElasticAnisotropic::DeliverClone( string const & name,
                                      ManagedGroup * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<LinearElasticAnisotropic>( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  LinearElasticAnisotropic * const newConstitutiveRelation = dynamic_cast<LinearElasticAnisotropic *>(clone.get());


  newConstitutiveRelation->m_defaultStiffness = m_defaultStiffness;
  newConstitutiveRelation->m_stiffness = m_stiffness;
}

void LinearElasticAnisotropic::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_stiffness.resize( parent->size() );
  m_stiffness = m_defaultStiffness;
}

void LinearElasticAnisotropic::PostProcessInput()
{
  m_stiffness = m_defaultStiffness;
}

void LinearElasticAnisotropic::StateUpdatePoint( localIndex const k,
                                                 localIndex const q,
                                                 R2SymTensor const & D,
                                                 R2Tensor const & Rot,
                                                 integer const updateStiffnessFlag )
{
  R2SymTensor T;
  real64 * const restrict Tdata = T.Data();
  real64 const * const restrict Ddata = D.Data();
  real64 const (&c)[6][6] = m_stiffness[k].m_data;

  Tdata[0] += c[0][0]*Ddata[0] + c[0][1]*Ddata[2] + c[0][2]*Ddata[5] + c[0][3]*(2*Ddata[4]) + c[0][4]*(2*Ddata[3]) + c[0][5]*(2*Ddata[1]);
  Tdata[2] += c[1][0]*Ddata[0] + c[1][1]*Ddata[2] + c[1][2]*Ddata[5] + c[1][3]*(2*Ddata[4]) + c[1][4]*(2*Ddata[3]) + c[1][5]*(2*Ddata[1]);
  Tdata[5] += c[2][0]*Ddata[0] + c[2][1]*Ddata[2] + c[2][2]*Ddata[5] + c[2][3]*(2*Ddata[4]) + c[2][4]*(2*Ddata[3]) + c[2][5]*(2*Ddata[1]);
  Tdata[4] += c[3][0]*Ddata[0] + c[3][1]*Ddata[2] + c[3][2]*Ddata[5] + c[3][3]*(2*Ddata[4]) + c[3][4]*(2*Ddata[3]) + c[3][5]*(2*Ddata[1]);
  Tdata[3] += c[4][0]*Ddata[0] + c[4][1]*Ddata[2] + c[4][2]*Ddata[5] + c[4][3]*(2*Ddata[4]) + c[4][4]*(2*Ddata[3]) + c[4][5]*(2*Ddata[1]);
  Tdata[1] += c[5][0]*Ddata[0] + c[5][1]*Ddata[2] + c[5][2]*Ddata[5] + c[5][3]*(2*Ddata[4]) + c[5][4]*(2*Ddata[3]) + c[5][5]*(2*Ddata[1]);


  m_meanStress[k][q] = ( Tdata[0] + Tdata[2] + Tdata[5] ) / 3.0;
  T.PlusIdentity( -m_meanStress[k][q] );
  m_deviatorStress[k][q] += T;

  T.QijAjkQlk( m_deviatorStress[k][q], Rot );
  m_deviatorStress[k][q] = T;

}

void LinearElasticAnisotropic::GetStiffness( localIndex const k, real64 c[6][6] ) const
{
  for( int i=0 ; i<6 ; ++i )
  {
    for( int j=0 ; j<6 ; ++j )
    {
      c[i][j] = m_stiffness[k].m_data[i][j];
    }
  }
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticAnisotropic, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
