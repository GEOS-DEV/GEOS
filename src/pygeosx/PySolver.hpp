#include "events/EventBase.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "mainInterface/ProblemManager.hpp"


namespace geosx
{
namespace python
{


template<typename T>
PyObject * explicitStep(T * self, PyObject * args)
{
  //VERIFY_NON_NULL_SELF( self );
  //VERIFY_INITIALIZED( self );

  double time;
  double dt;
  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geosx::DomainPartition & domain = self->pb_manager->getDomainPartition();

  self->group->explicitStep(time, dt, 0, domain);

  Py_RETURN_NONE;
}

PyObject * createNewPySolver( geosx::EventBase * subEvent, geosx::ProblemManager * group );

PyTypeObject * getPySolverType();

} // namespace python
} // namespace geosx
