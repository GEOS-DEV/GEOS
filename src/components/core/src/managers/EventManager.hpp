/*
 * EventManager.hpp
 *
 *  Created on: Oct 5, 2016
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"


namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const Events("Events");
}
}

class EventManager : public dataRepository::ManagedGroup
{
public:
  EventManager( std::string const & name,
                ManagedGroup * const parent );

  virtual ~EventManager() override;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  void CheckEventTiming();

};

class SolverApplication : public dataRepository::ManagedGroup
{
public:
  SolverApplication( std::string const & name,
                     ManagedGroup * const parent );

  virtual ~SolverApplication() override;

  static string CatalogName() { return "EventManager"; }

  virtual void FillDocumentationNode() override;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_ */
