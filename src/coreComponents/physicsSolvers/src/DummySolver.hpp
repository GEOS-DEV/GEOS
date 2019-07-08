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


#ifndef DUMMYSOLVER_HPP_
#define DUMMYSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class DomainPartition;

class DummySolver : public SolverBase
{
public:
  DummySolver( const std::string& name,
                           ManagedGroup * const parent );


  virtual ~DummySolver() override;

  static string CatalogName() { return "DummySolver"; }


  virtual real64 SolverStep( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         DomainPartition * domain ) override;

  virtual real64 GetTimestepRequest(real64 const time) override;

  struct viewKeyStruct
  {
    static constexpr auto randScaleString = "scale";
    static constexpr auto randSeedString = "seed";
    dataRepository::ViewKey randScale = { "scale" };
    dataRepository::ViewKey randSeed = { "seed" };
    } viewKeysDummySolver;

protected:
  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override final;

  real64 m_randScale;
  integer m_randSeed;

};

} /* namespace geosx */

#endif /* DUMMYSOLVER_HPP_ */
