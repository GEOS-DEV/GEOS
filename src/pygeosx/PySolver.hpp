#include "events/EventBase.hpp"
#include "physicsSolvers/SolverBase.hpp"


namespace geosx
{
namespace python
{

struct PySolver
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geosx::SolverBase.";

  geosx::SolverBase * solver;
};


PyObject * explicitStep(PySolver * self, PyObject * args);

PyObject * createNewPySolver( geosx::EventBase * subEvent );

PyTypeObject * getPySolverType();

} // namespace python
} // namespace geosx
