#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "pygeosx.hpp"
#include "PySolver.hpp"

// SolverBase

#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PySolver is not initialized.", nullptr )


namespace geosx
{

namespace python
{

struct PySolver
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geosx::SolverBase.";

  geosx::SolverBase * group;
  geosx::ProblemManager * pb_manager;
};


static PyObject * PySolver_repr( PyObject * const obj ) noexcept
{
  PySolver const * const pySolver = LvArray::python::convert< PySolver >( obj, getPySolverType() );
  if( pySolver == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pySolver );

  string const path = pySolver->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pySolver->group) ).name() );
  string const repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );
}




BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef PySolver_methods[] = {
{ "explicitStep", (PyCFunction) explicitStep<PySolver>, METH_VARARGS, "explicit Step" },
{ nullptr, nullptr, 0, nullptr }        /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */
static PyTypeObject PySolverType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.Solver",
  .tp_basicsize = sizeof( PySolver ),
  .tp_itemsize = 0,
  .tp_repr = PySolver_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = PySolver::docString,
  .tp_methods = PySolver_methods,
  .tp_new = PyType_GenericNew,
};


END_ALLOW_DESIGNATED_INITIALIZERS


PyObject * createNewPySolver( geosx::EventBase * subEvent, geosx::ProblemManager * pbManager )
{
  // Create a new Group and set the dataRepository::Group it points to.
  PyObject * const ret = PyObject_CallFunction( reinterpret_cast< PyObject * >( getPySolverType() ), "" );
  PySolver * const retSolver = reinterpret_cast< PySolver * >( ret );
  if( retSolver == nullptr )
  {
    return nullptr;
  }

  geosx::SolverBase * solver = static_cast<geosx::SolverBase *>(subEvent->getEventTarget());
  retSolver->group = solver;
  retSolver->pb_manager = pbManager;

  return ret;
}


PyTypeObject * getPySolverType()
{ return &PySolverType; }

}
}
