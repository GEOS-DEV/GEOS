/*
 * SolverBase.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#include "SolverBase.hpp"

namespace geosx
{

SolverBase::SolverBase( std::string const & name,
                        WrapperCollection * const parent ) :
  WrapperCollection( name, parent )
{
  // TODO Auto-generated constructor stub

}

SolverBase::~SolverBase()
{
  // TODO Auto-generated destructor stub
}

} /* namespace ANST */
