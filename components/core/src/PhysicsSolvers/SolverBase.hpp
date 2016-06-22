/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_

#include "DataTypes.hpp"
#include "ObjectCatalogue.hpp"
#include <string>

namespace geosx
{
class DataObjectManager;

class SolverBase
{
public:
  typedef ObjectCatalogueEntryBase< geosx::SolverBase, std::string > SolverFactory;

  SolverBase( std::string const & name );

  virtual ~SolverBase();


  virtual void RegisterDataObjects( DataObjectManager& domain ) = 0;

  virtual void TimeStep( const real64& time_n,
                         const real64& dt,
                         const int cycleNumber,
                         DataObjectManager& domain ) = 0;


/*
  static std::map< std::string, SolverBase* >& GetCatalogue()
  {
    static std::map< std::string, SolverBase* > m_catalogue;
    return m_catalogue;
  }

  template< typename DERIVED >
  void addCatalogueEntry()
  {
    std::string name = DERIVED::CatalogueName();
    SolverBase::GetCatalogue()[name] = this;
  }

  template< typename DERIVED, typename ...ARGS >
  static std::unique_ptr<SolverBase> Factory( const std::string& name, const ARGS&... args )
  {
    return std::unique_ptr<SolverBase>( new DERIVED(args...) );
  }

*/

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
