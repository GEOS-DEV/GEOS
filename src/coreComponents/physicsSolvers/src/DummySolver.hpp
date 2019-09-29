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


#ifndef DUMMYSOLVER_HPP_
#define DUMMYSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
namespace dataRepository
{
class Group;
}
class DomainPartition;

class DummySolver : public SolverBase
{
public:
  DummySolver( const std::string& name,
                           Group * const parent );


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
  virtual void InitializePreSubGroups( Group * const problemManager ) override final;

  real64 m_randScale;
  integer m_randSeed;

};

} /* namespace geosx */

#endif /* DUMMYSOLVER_HPP_ */
