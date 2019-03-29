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
  m_stiffness0{},
  m_stiffness{}
{
  RegisterViewWrapper( viewKeyStruct::stiffness0String, &m_stiffness0, 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Elastic Stiffness Tensor");

  RegisterViewWrapper( viewKeyStruct::stiffnessString, &m_stiffness, 0 )->
    setApplyDefaultValue(-1)->
    setDescription("Elastic Modulus Field");

}


LinearElasticAnisotropic::~LinearElasticAnisotropic()
{}


void
LinearElasticAnisotropic::DeliverClone( string const & name,
                                      ManagedGroup * const parent,
                                      std::unique_ptr<ConstitutiveBase> & clone ) const
{
  std::unique_ptr<LinearElasticAnisotropic>
  newConstitutiveRelation = std::make_unique<LinearElasticAnisotropic>( name, parent );

  newConstitutiveRelation->m_meanStress = m_meanStress;
  newConstitutiveRelation->m_deviatorStress = m_deviatorStress;

  clone = std::move(newConstitutiveRelation);
}

void LinearElasticAnisotropic::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_stiffness.resize( parent->size() );
}

void LinearElasticAnisotropic::PostProcessInput()
{

}

void LinearElasticAnisotropic::StateUpdatePoint( localIndex const i,
                                               localIndex const q,
                                               R2SymTensor const & D,
                                               R2Tensor const & Rot,
                                               integer const systemAssembleFlag )
{

}

void LinearElasticAnisotropic::GetStiffness( localIndex const k, real64 c[6][6] ) const
{
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearElasticAnisotropic, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
