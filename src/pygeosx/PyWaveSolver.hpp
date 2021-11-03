#include "events/EventBase.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{
namespace python
{

PyObject * createNewPyWaveSolver( geosx::EventBase * subEvent, geosx::ProblemManager * pbManager );

PyTypeObject * getPyWaveSolverType();

} // namespace python
} // namespace geosx
