/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file ConstitutiveBase.hpp
 */

#ifndef CONSTITUTIVEBASE_HPP_
#define CONSTITUTIVEBASE_HPP_

#include "common/DataTypes.hpp"
#include "ObjectCatalog.hpp"
#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const stateData( "StateData" );
string const parameterData( "ParameterData" );
}
}

namespace constitutive
{


class ConstitutiveBase : public dataRepository::ManagedGroup
{
public:


  ConstitutiveBase( std::string const & name,
                    ManagedGroup * const parent );

  virtual ~ConstitutiveBase() override;

  virtual std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                          ManagedGroup * const parent ) const = 0;


  virtual void SetParamStatePointers( void *& ) = 0;


  typedef void (*UpdateFunctionPointer)( R2SymTensor const & D,
                                         R2Tensor const & Rot,
                                         localIndex const i,
                                         void * dataPtrs,
                                         integer const systemAssembleFlag );

  virtual UpdateFunctionPointer GetStateUpdateFunctionPointer( ) = 0;

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const = 0;


  virtual R2SymTensor StateUpdatePoint( R2SymTensor const & D,
                                        R2Tensor const & Rot,
                                        localIndex const i,
                                        localIndex const q,
                                        integer const systemAssembleFlag ) { return R2SymTensor(); }

  virtual void FluidPressureUpdate( real64 const &dens,
                                    localIndex const i,
                                    real64 &pres,
                                    real64 &dPres_dDens ) {}

  virtual void FluidDensityUpdate( real64 const &pres,
                                   localIndex const i,
                                   real64 &dens,
                                   real64 &dDens_dPres ) {}

  virtual void DensityUpdate( real64 const &pres,
                              localIndex const k,
                              localIndex const q ) {}


  virtual void FluidViscosityUpdate( real64 const &pres,
                                     localIndex const i,
                                     real64 &visc,
                                     real64 &dVisc_dPres ) {}

  virtual void SimplePorosityUpdate( real64 const &pres,
                                     real64 const &poro_ref,
                                     localIndex const i,
                                     real64 &poro,
                                     real64 &dPoro_dPres ) {}

  virtual void FillDocumentationNode() override = 0;

  virtual void resize( localIndex ) override;

  virtual void GetStiffness( realT c[6][6] ) const = 0;


  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, ManagedGroup * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

  virtual string GetCatalogName() = 0;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex );

  struct viewKeyStruct
  {} m_ConstitutiveBaseViewKeys;

  struct groupKeyStruct
  {} m_ConstitutiveBaseGroupKeys;

  virtual viewKeyStruct       & viewKeys()        { return m_ConstitutiveBaseViewKeys; }
  virtual viewKeyStruct const & viewKeys() const  { return m_ConstitutiveBaseViewKeys; }

  virtual groupKeyStruct       & groupKeys()       { return m_ConstitutiveBaseGroupKeys; }
  virtual groupKeyStruct const & groupKeys() const { return m_ConstitutiveBaseGroupKeys; }

protected:

private:
  ManagedGroup * m_constitutiveDataGroup = nullptr;

  ConstitutiveBase( ConstitutiveBase const & ) = delete;
  ConstitutiveBase( ConstitutiveBase && ) = delete;
  ConstitutiveBase const & operator=( ConstitutiveBase const & ) = delete;
  ConstitutiveBase const & operator=( ConstitutiveBase && ) = delete;

};



}
}
#endif
