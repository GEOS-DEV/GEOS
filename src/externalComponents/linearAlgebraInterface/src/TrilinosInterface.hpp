/**
 * @file TrilinosInterface.hpp
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "EpetraSparseMatrix.hpp"
#include "EpetraVector.hpp"
#include "TrilinosSolver.hpp"

namespace geosx
{

class TrilinosInterface
{
public:

  // Epetra matrix and vector wrappers
  using ParallelMatrix = EpetraSparseMatrix;
  using ParallelVector = EpetraVector;

  using LinearSolver = TrilinosSolver;

  TrilinosInterface() = default;
  ~TrilinosInterface() = default;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_ */
