#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "pygeosx.hpp"
#include "PySolver.hpp"
#include "PyWaveSolver.hpp"
#include "PyGroup.hpp"

#include "physicsSolvers/wavePropagation/WaveSolverBase.hpp"

#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PySolver is not initialized.", nullptr )

namespace geosx
{

namespace python
{


struct PyWaveSolver
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geosx::WaveSolverBase.";

  geosx::WaveSolverBase * group;
  geosx::ProblemManager * pb_manager;
};


static PyObject * PyWaveSolver_repr( PyObject * const obj ) noexcept
{
  PyWaveSolver const * const pyWaveSolver = LvArray::python::convert< PyWaveSolver >( obj, getPyWaveSolverType() );
  if( pyWaveSolver == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyWaveSolver );

  string const path = pyWaveSolver->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pyWaveSolver->group) ).name() );
  string const repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );
}



static PyObject * postProcessInput(PyWaveSolver * self, PyObject * args) noexcept
{
  GEOSX_UNUSED_VAR( self, args);

  self->group->postProcessInput();

  Py_RETURN_NONE;
}



static PyObject * precomputeSourceAndReceiverTerm(PyWaveSolver * self, PyObject * args) noexcept
{
  GEOSX_UNUSED_VAR( self, args);

  geosx::DomainPartition & domain = self->pb_manager->getDomainPartition();
  geosx::MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  self->group->precomputeSourceAndReceiverTerm(mesh);

  Py_RETURN_NONE;
}


BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef PyWaveSolver_methods[] = {
{ "explicitStep", (PyCFunction) explicitStep<PyWaveSolver>, METH_VARARGS, "explicit Step" },
{ "resizeArrays", (PyCFunction) postProcessInput, METH_NOARGS, "resize pressure_at_receivers array to fit with new number of receivers"},
{ "updateSourceAndReceivers", (PyCFunction) precomputeSourceAndReceiverTerm, METH_NOARGS, "update positions for new source and receivers"},
{ "get_wrapper", (PyCFunction) PyGroup_getWrapper<PyWaveSolver>, METH_VARARGS, PyGroup_getWrapperDocString },
{ nullptr, nullptr, 0, nullptr }        /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */
static PyTypeObject PyWaveSolverType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.WaveSolver",
  .tp_basicsize = sizeof( PyWaveSolver ),
  .tp_itemsize = 0,
  .tp_repr = PyWaveSolver_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = PyWaveSolver::docString,
  .tp_methods = PyWaveSolver_methods,
  .tp_new = PyType_GenericNew,
};


END_ALLOW_DESIGNATED_INITIALIZERS


PyObject * createNewPyWaveSolver( geosx::EventBase * subEvent, geosx::ProblemManager * pbManager )
{
  PyObject * ret = PyObject_CallFunction(reinterpret_cast< PyObject * >( getPyWaveSolverType() ), "");
  PyWaveSolver * const retSolver = reinterpret_cast< PyWaveSolver * >( ret );

  if( retSolver == nullptr )
  {
    return nullptr;
  }

  geosx::WaveSolverBase * solver = static_cast<geosx::WaveSolverBase *>(subEvent->getEventTarget());
  retSolver->group = solver;
  retSolver->pb_manager = pbManager;


  return ret;
}


PyTypeObject * getPyWaveSolverType()
{ return &PyWaveSolverType; }

}
}
