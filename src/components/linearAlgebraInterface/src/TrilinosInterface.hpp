/**
 * @file TrilinosInterface.hpp
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "EpetraSparseMatrix.hpp"
#include "EpetraVector.hpp"

namespace geosx
{


class TrilinosInterface
{
public:

  // your wrappers go here instead of the naked epetra types
  using ParallelMap = Epetra_Map;
  using ParallelGraph = Epetra_CrsGraph;
  using ParallelMatrix = EpetraSparseMatrix;
  using ParallelVector = EpetraVector;


  TrilinosInterface() = default;
  ~TrilinosInterface() = default;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_ */
