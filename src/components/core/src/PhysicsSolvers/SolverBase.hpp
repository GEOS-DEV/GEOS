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

#include "../dataRepository/SynchronizedGroup.hpp"
#include "common/DataTypes.hpp"
#include <pugixml.hpp>


namespace geosx
{

class SolverBase : public dataRepository::SynchronizedGroup
{
public:

  explicit SolverBase( std::string const & name,
                       SynchronizedGroup * const parent );

  virtual ~SolverBase();

  SolverBase() = default;
  SolverBase( SolverBase const & ) = default;
  SolverBase( SolverBase &&) = default;
  SolverBase& operator=( SolverBase const & ) = default;
  SolverBase& operator=( SolverBase&& ) = default;


//  virtual void Registration( dataRepository::WrapperCollection& domain );

  virtual void ReadXML( pugi::xml_node solverNode );

  virtual void Registration( dataRepository::SynchronizedGroup * const domain ) override;

  virtual void Initialize( dataRepository::SynchronizedGroup& domain );

  virtual void TimeStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::SynchronizedGroup& domain ) = 0;



  using CatalogInterface = cxx_utilities::CatalogInterface< SolverBase, std::string const &, SynchronizedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

private:

};



} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
