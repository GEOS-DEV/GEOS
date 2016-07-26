/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

#include "dataRepository/intrinsic/WrapperCollection.hpp"

namespace geosx
{

class ProblemManager : public dataRepository::WrapperCollection
{
public:
  explicit ProblemManager();
  ~ProblemManager();

  void ParseCommandLineInput( int const& argc, char* const argv[]) ;

  void ParseInputFile();

};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
