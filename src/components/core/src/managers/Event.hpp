/*
 * Event.hpp
 *
 *  Created on: Jul 22, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENT_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENT_HPP_

namespace geosx
{

class Event : public dataRepository::ManagedGroup
{
public:
  Event() = delete;
  Event( string & name, ManagedGroup * parent );
  virtual ~Event();

  virtual void Execute(  ) = 0;



};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENT_HPP_ */
