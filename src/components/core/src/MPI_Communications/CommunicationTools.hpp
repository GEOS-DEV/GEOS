/*
 * CommunicationTools.hpp
 *
 *  Created on: Jan 6, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_
#define SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_

#include "common/DataTypes.hpp"
#include "mpi.h"
namespace geosx
{


class ObjectManagerBase;
class NeighborCommunicator;

class CommunicationTools
{
public:
  CommunicationTools();
  ~CommunicationTools();

  static void AssignGlobalIndices( ObjectManagerBase & object,
                                   ObjectManagerBase const & compositionObject,
                                   array<NeighborCommunicator> & neighbors );

  static int MPI_Size( MPI_Comm const & comm );
  static int MPI_Rank( MPI_Comm const & comm );
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_ */
