#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "physicsSolvers/SolverBase.hpp"
#include "LvArray/src/python/pythonForwardDeclarations.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"

#include "PySolver.hpp"
#include "dataRepository/python/PyGroup.hpp"

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
};


static PyObject * PySolver_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOSX_UNUSED_VAR( args, kwds );
  PySolver *self;

  self = (PySolver *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->group = nullptr;
  }

  return (PyObject *)self;
}


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



static PyObject * solverStep( PySolver * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;
  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geosx::GeosxState * g_state = &getGlobalState();
  geosx::DomainPartition & domain = g_state->getProblemManager().getDomainPartition();

  self->group->solverStep( time, dt, 0, domain );

  Py_RETURN_NONE;
}


static PyObject * reinit( PySolver * self, PyObject *args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOSX_UNUSED_VAR( args );

  self->group->reinit();

  Py_RETURN_NONE;
}


static PyMethodDef PySolver_methods[] = {
  { "solverStep", (PyCFunction) solverStep, METH_VARARGS, "solver Step" },
  { "reinit", (PyCFunction) reinit, METH_NOARGS, "re-initialize certain variable depending on the solver being used"},
  { "get_wrapper", (PyCFunction) PyGroup_getWrapper< PySolver >, METH_VARARGS, PyGroup_getWrapperDocString },
  { nullptr, nullptr, 0, nullptr }      /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PySolverType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pysolver.Solver",
  .tp_basicsize = sizeof( PySolver ),
  .tp_itemsize = 0,
  .tp_repr = PySolver_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PySolver::docString,
  .tp_methods = PySolver_methods,
  .tp_new = PySolver_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef SolverFuncs[] = {
  { nullptr, nullptr, 0, nullptr }        /* Sentinel */
};

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyModuleDef pysolvermodule = {
  PyModuleDef_HEAD_INIT,
  .m_name = "pysolver",
  .m_doc = "pysolver module for SolverBase",
  .m_size = -1,
  .m_methods = SolverFuncs
};

END_ALLOW_DESIGNATED_INITIALIZERS

PyMODINIT_FUNC
PyInit_pysolver( void )
{
  LvArray::python::PyObjectRef<> module{ PyModule_Create( &pysolvermodule ) };
  if( module == nullptr )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geosx::python::getPySolverType(), "Solver" ) )
  {
    return nullptr;
  }

  return module.release();
}

PyTypeObject * getPySolverType()
{ return &PySolverType; }

}
}
