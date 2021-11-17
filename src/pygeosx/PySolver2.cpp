#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "pygeosx.hpp"
#include "PySolver2.hpp"
#include "PyGroup.hpp"

#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PySolver is not initialized.", nullptr )


namespace geosx
{

namespace python
{

struct PySolver2
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geosx::SolverBase.";

  geosx::SolverBase * group;
  geosx::ProblemManager * pb_manager;
};


static PyObject * PySolver2_repr( PyObject * const obj ) noexcept
{
  PySolver2 const * const pySolver = LvArray::python::convert< PySolver2 >( obj, getPySolver2Type() );
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



static PyObject * postProcessInput(PySolver2* self, PyObject* args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOSX_UNUSED_VAR( args);

  self->group->postProcessInput();

  Py_RETURN_NONE;
}


static PyObject * initPostInitialConditions(PySolver2* self, PyObject* args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOSX_UNUSED_VAR( args);

  self->group->initializePostInitialConditionsPreSubGroups();

  Py_RETURN_NONE;
}


BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyMethodDef PySolver2_methods[] = {
{ "explicitStep", (PyCFunction) explicitStep<PySolver2>, METH_VARARGS, "explicit Step" },
{ "linearImplicitStep", (PyCFunction) linearImplicitStep<PySolver2>, METH_VARARGS, "linear implicit step" },
{ "nonlinearImplicitStep", (PyCFunction) nonlinearImplicitStep<PySolver2>, METH_VARARGS, "non linear implicit step" },
{ "postProcessInput", (PyCFunction) postProcessInput, METH_NOARGS, "post processing input"},
{ "initPostInitialConditions", (PyCFunction) initPostInitialConditions, METH_NOARGS, "call initializePostInitialConditionsPreSubGroup"},
{ "get_wrapper", (PyCFunction) PyGroup_getWrapper<PySolver2>, METH_VARARGS, PyGroup_getWrapperDocString },
{ nullptr, nullptr, 0, nullptr }        /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */
static PyTypeObject PySolver2Type = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.Solver",
  .tp_basicsize = sizeof( PySolver2 ),
  .tp_itemsize = 0,
  .tp_repr = PySolver2_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = PySolver2::docString,
  .tp_methods = PySolver2_methods,
  .tp_new = PyType_GenericNew,
};


END_ALLOW_DESIGNATED_INITIALIZERS


PyObject * createNewPySolver2( geosx::EventBase * subEvent, geosx::ProblemManager * pbManager )
{
  // Create a new Group and set the dataRepository::Group it points to.
  PyObject * const ret = PyObject_CallFunction( reinterpret_cast< PyObject * >( getPySolver2Type() ), "" );
  PySolver2 * const retSolver = reinterpret_cast< PySolver2 * >( ret );
  if( retSolver == nullptr )
  {
    return nullptr;
  }

  geosx::SolverBase * solver = static_cast<geosx::SolverBase *>(subEvent->getEventTarget());
  retSolver->group = solver;
  retSolver->pb_manager = pbManager;

  return ret;
}


PyTypeObject * getPySolver2Type()
{ return &PySolver2Type; }

}
}
