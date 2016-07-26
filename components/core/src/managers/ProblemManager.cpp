/*
 * ProblemManager.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#include "ProblemManager.hpp"

namespace geosx
{

ProblemManager::ProblemManager():
    WrapperCollection("ProblemManager",nullptr)
{}

ProblemManager::~ProblemManager()
{}

void ProblemManager::ParseCommandLineInput( int const& argc, char* const argv[])
{

}


} /* namespace geosx */
