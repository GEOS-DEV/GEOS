#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "pygeosx.hpp"
#include "PySolver.hpp"

// SolverBase
#include "physicsSolvers/WaveSolverBase.hpp"

#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->solver == nullptr, PyExc_RuntimeError, "The PySolver is not initialized.", nullptr )

namespace geosx
{

namespace python
{


struct PyWaveSolver
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geosx::WaveSolverBase.";

  geosx::WaveSolverBase * solver;
};


static PyObject * PyWaveSolver_repr( PyObject * const obj ) noexcept
{
  PySolver const * const pyWaveSolver = LvArray::python::convert< PyWaveSolver >( obj, getPySolverType() );
  if( pySolver == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyWaveSolver );

  string const path = pyWaveSolver->solver->getPath();
  string const type = LvArray::system::demangle( typeid( *(pyWaveSolver->solver) ).name() );
  string const repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );
}



static PyObject * postProcessInput(PyWaveSolver * self, PyObject * args) noexcept
{
  GEOSX_UNUSED_VAR( self, args);

  self->solver->postProcessInput();

  Py_RETURN_NONE;
}



static PyObject * precomputeSourceAndReceiverTerm(PyWaveSolver * self, PyObject * args) noexcept
{
  GEOSX_UNUSED_VAR( self, args);

  geosx::DomainPartition & domain = geosx::g_state->getProblemManager().getDomainPartition();
  geosx::MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  self->solver->precomputeSourceAndReceiverTerm(mesh);

  Py_RETURN_NONE;
}


BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef PyWaveSolver_methods[] = {
{ "explicitStep", (PyCFunction) explicitStep, METH_VARARGS, "explicit Step" },
{ "postProcessInput", (PyCFunction) postProcessInput, METH_NOARGS, "resize pressure_at_receivers array to fit with new number of receivers"},
{ "precomputeSourceAndReceiverTerm", (PyCFunction) precomputeSourceAndReceiverTerm, METH_NOARGS, "update positions for new source and receivers"},
{ nullptr, nullptr, 0, nullptr }        /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */
static PyTypeObject PyWaveSolverType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.WaveSolvers",
  .tp_basicsize = sizeof( PyWaveSolver ),
  .tp_itemsize = 0,
  .tp_repr = PyWaveSolver_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = PyWaveSolver::docString,
  .tp_methods = PyWaveSolver_methods,
  .tp_new = PyType_GenericNew,
};


END_ALLOW_DESIGNATED_INITIALIZERS


PyObject * createNewPyWaveSolver( geosx::EventBase * subEvent )
{
  // Create a new Group and set the dataRepository::Group it points to.
  PyObject * const ret = PyObject_CallFunction( reinterpret_cast< PyObject * >( getPyWaveSolverType() ), "" );
  PySolver * const retSolver = reinterpret_cast< PySolver * >( ret );
  if( retSolver == nullptr )
  {
    return nullptr;
  }

  geosx::WaveSolverBase * solver = static_cast<geosx::WaveSolverBase *>(subEvent->getEventTarget());
  retSolver->solver = solver;

  return ret;
}


PyTypeObject * getPyWaveSolverType()
{ return &PyWaveSolverType; }

}
}
