/**
 * @file TrilinosInterface.hpp
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

namespace geosx
{


class TrilinosInterface
{
public:

  // your wrappers go here instead of the naked epetra types
  using ParallelMap = Epetra_Map;
  using ParallelGraph = Epetra_FECrsGraph;
  using ParallelMatrix = Epetra_FECrsMatrix;
  using ParallelVector = Epetra_FEVector;


  TrilinosInterface() = default;
  ~TrilinosInterface() = default;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_ */
