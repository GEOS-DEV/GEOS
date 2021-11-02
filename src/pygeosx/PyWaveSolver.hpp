#include "events/EventBase.hpp"

namespace geosx
{
namespace python
{

PyObject * createNewPySolver( geosx::EventBase * subEvent );

PyTypeObject * getPySolverType();

} // namespace python
} // namespace geosx
