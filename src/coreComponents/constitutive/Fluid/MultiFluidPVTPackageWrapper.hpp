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
  * @file MultiFluidPVTPackageWrapper.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDPVTPACKAGEWRAPPER_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDPVTPACKAGEWRAPPER_HPP

#include "constitutive/Fluid/MultiFluidBase.hpp"

// PVTPackage includes
#include "MultiphaseSystem/PVTEnums.hpp"

namespace PVTPackage
{
class MultiphaseSystem;
}

namespace geosx
{

namespace constitutive
{

class MultiFluidPVTPackageWrapper : public MultiFluidBase
{
public:

  MultiFluidPVTPackageWrapper( std::string const & name, ManagedGroup * const parent );

  virtual ~MultiFluidPVTPackageWrapper() override;

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override final {}



  virtual void StateUpdatePointMultiFluid( real64 const & pres,
                                           real64 const & temp,
                                           arraySlice1d<real64 const> const & composition,
                                           localIndex const k,
                                           localIndex const q ) override;



protected:
  virtual void PostProcessInput() override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;

  // function that populates m_fluid ptr; to be overriden by derived classes
  virtual void createFluid() = 0;

  // PVTPackage phase type labels
  array1d<PVTPackage::PHASE_TYPE> m_pvtPackagePhaseTypes;

  // PVTPackage fluid object
  PVTPackage::MultiphaseSystem * m_fluid;

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDPVTPACKAGEWRAPPER_HPP
