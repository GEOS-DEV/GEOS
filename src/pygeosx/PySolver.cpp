#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "pygeosx.hpp"
#include "PySolver.hpp"

// SolverBase

#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->solver == nullptr, PyExc_RuntimeError, "The PySolver is not initialized.", nullptr )

namespace geosx
{

namespace python
{



static PyObject * PySolver_repr( PyObject * const obj ) noexcept
{
  PySolver const * const pySolver = LvArray::python::convert< PySolver >( obj, getPySolverType() );
  if( pySolver == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pySolver );

  string const path = pySolver->solver->getPath();
  string const type = LvArray::system::demangle( typeid( *(pySolver->solver) ).name() );
  string const repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );
}



PyObject * explicitStep(PySolver * self, PyObject * args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;
  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geosx::DomainPartition & domain = geosx::g_state->getProblemManager().getDomainPartition();

  self->solver->explicitStep(time, dt, 0, domain);

  Py_RETURN_NONE;
}



BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef PySolver_methods[] = {
{ "explicitStep", (PyCFunction) explicitStep, METH_VARARGS, "explicit Step" },
{ nullptr, nullptr, 0, nullptr }        /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */
static PyTypeObject PySolverType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.Solvers",
  .tp_basicsize = sizeof( PySolver ),
  .tp_itemsize = 0,
  .tp_repr = PySolver_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = PySolver::docString,
  .tp_methods = PySolver_methods,
  .tp_new = PyType_GenericNew,
};


END_ALLOW_DESIGNATED_INITIALIZERS


PyObject * createNewPySolver( geosx::EventBase * subEvent )
{
  // Create a new Group and set the dataRepository::Group it points to.
  PyObject * const ret = PyObject_CallFunction( reinterpret_cast< PyObject * >( getPySolverType() ), "" );
  PySolver * const retSolver = reinterpret_cast< PySolver * >( ret );
  if( retSolver == nullptr )
  {
    return nullptr;
  }

  geosx::SolverBase * solver = static_cast<geosx::SolverBase *>(subEvent->getEventTarget());
  retSolver->solver = solver;

  return ret;
}


PyTypeObject * getPySolverType()
{ return &PySolverType; }

}
}
