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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_



#include <string>
#include <limits>

#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "../dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"
#include "mesh/MeshBody.hpp"


namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const courant = "courant";
string const maxDt   = "maxDt";
}
}

class SolverBase : public dataRepository::ManagedGroup
{
public:

  explicit SolverBase( std::string const & name,
                       ManagedGroup * const parent );

  virtual ~SolverBase() override;

  static string CatalogName() { return "SolverBase"; }

  SolverBase() = default;
  SolverBase( SolverBase const & ) = default;
  SolverBase( SolverBase &&) = default;
  SolverBase& operator=( SolverBase const & ) = default;
  SolverBase& operator=( SolverBase&& ) = default;


//  virtual void Registration( dataRepository::WrapperCollection& domain );


  virtual void TimeStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::ManagedGroup * domain );

  // virtual void TimeStep( real64 const & time_n,
  //                        real64 const & dt,
  //                        int const cycleNumber,
  //                        dataRepository::ManagedGroup * domain ) = 0;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  using CatalogInterface = cxx_utilities::CatalogInterface< SolverBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  struct viewKeyStruct
  {
    constexpr static auto verboseLevelString = "verboseLevel";
    constexpr static auto gravityVectorString = "gravityVector";

  } viewKeys;

  struct groupKeyStruct
  {
  } groupKeys;

  R1Tensor const & getGravityVector() const { return m_gravityVector; }
  R1Tensor       & getGravityVector()       { return m_gravityVector; }
  R1Tensor const * globalGravityVector() const;
  integer verboseLevel() const { return m_verboseLevel; }
private:
  integer m_verboseLevel = 0;
  R1Tensor m_gravityVector;

};



} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
