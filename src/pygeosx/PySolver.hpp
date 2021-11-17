#include "events/EventBase.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "mainInterface/ProblemManager.hpp"


namespace geosx
{
namespace python
{


PyMODINIT_FUNC PyInit_pysolver(void);

PyTypeObject * getPySolverType();



} // namespace python
} // namespace geosx
