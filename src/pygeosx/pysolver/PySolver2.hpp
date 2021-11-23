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



template<typename T>
PyObject * linearImplicitStep(T * self, PyObject * args)
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

  self->group->linearImplicitStep(time, dt, 0, domain);

  Py_RETURN_NONE;
}



template<typename T>
PyObject * nonlinearImplicitStep(T * self, PyObject * args)
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

  self->group->nonlinearImplicitStep(time, dt, 0, domain);

  Py_RETURN_NONE;
}

PyObject * createNewPySolver2( geosx::EventBase * subEvent, geosx::ProblemManager * group );

PyTypeObject * getPySolver2Type();



} // namespace python
} // namespace geosx
