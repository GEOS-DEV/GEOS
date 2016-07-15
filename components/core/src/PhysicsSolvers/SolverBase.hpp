/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_

#include "dataRepository/DataTypes.hpp"
#include "codingUtilities/ObjectCatalogue.hpp"
#include <string>

namespace geosx
{
namespace dataRepository
{
class WrapperCollection;
}

class SolverBase
{
public:
//  typedef ObjectCatalogueEntryBase< geosx::SolverBase, std::string > SolverFactory;

  SolverBase( std::string const & name );

  virtual ~SolverBase();

  virtual void Registration( dataRepository::WrapperCollection& domain ) = 0;

  virtual void TimeStep( const real64& time_n,
                         const real64& dt,
                         const int cycleNumber,
                         dataRepository::WrapperCollection& domain ) = 0;

  CATALOGUE( SolverBase, VA_LIST( std::string const & name ), VA_LIST( name ) )

private:
  std::string m_name;

  SolverBase() = delete;
  SolverBase(const SolverBase&) = delete;
  SolverBase(const SolverBase&&) = delete;
  SolverBase& operator=(const SolverBase&) = delete;
  SolverBase& operator=(const SolverBase&&) = delete;
};





} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
