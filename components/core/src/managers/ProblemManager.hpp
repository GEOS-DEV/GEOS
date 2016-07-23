/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

namespace geosx
{

class ProblemManager : public dataRepository::WrapperCollection
{
public:
  explicit ProblemManager();
  ~ProblemManager();

  void ParseCommandLineInput( const int& argc, char* const argv[]) ;

  void ParseInputFile();

};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
