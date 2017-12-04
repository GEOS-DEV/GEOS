/*
 * Event.cpp
 *
 *  Created on: Jul 22, 2017
 *      Author: settgast
 */

#include "Event.hpp"

namespace geosx
{

Event::Event( string & name, ManagedGroup * parent ):
  ManagedGroup(name,parent)
{}

Event::~Event()
{}

} /* namespace geosx */
