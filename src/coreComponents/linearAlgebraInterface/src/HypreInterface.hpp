/**
 * @file HypreInterface.hpp
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_

//#include "HypreSparseMatrix.hpp"
//#include "HypreVector.hpp"
//#include "HypreSolver.hpp"

namespace geosx
{

class HypreInterface
{
public:

  // Epetra matrix and vector wrappers
  using ParallelMatrix = HypreSparseMatrix;
  using ParallelVector = HypreVector;

  using LinearSolver = HypreSolver;

  HypreInterface() = default;
  ~HypreInterface() = default;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_ */
