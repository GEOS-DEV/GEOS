.. _ObjectCatalog:


ObjectCatalog
=============

The "ObjectCatalog" is a collection of classes that acts as a statically initialized factory. 
It functions in a similar manner to a classic 
`factory method <https://en.wikibooks.org/wiki/C%2B%2B_Programming/Code/Design_Patterns#Abstract_Factory>`_, 
except that there is no maintained list of derived objects that is required to create new objects. 
In other words, there is no case-switch/if-elseif block to maintain.
Instead, the ``ObjectCatalog`` creates a "catalog" of derived objects using a ``std::unordered_map``.
The "catalog" is filled when new types are declared through the declaration of a helper class named
``CatalogEntryConstructor``.

The key functional  features of the "ObjectCatalog" concept may be summarized as:

* Anonymous addition of new objects to the catalog.
  Because we use a statically initialized singleton map object to store the catalog,
  no knowledge of the contents of the catalog is required in the main code.
  Therefore if a proprietary/sensitive catalog entry is desired, it is only required
  that the object definition be outside of the main repository and tied into the build
  system through some non-specific mechanism (i.e. a link in the src/externalComponents
  directory) and the catalog entry will be registered in the catalog without sharing 
  any knowledge of its existence.
  Then a proprietary input file may refer to the object to call for its creation.
* Zero maintenance catalog. 
  Again, because we use a singleton map to store the catalog, there is no updating of 
  code required to add new entries into the catalog. 
  The only modifications required are the actual source files of the catalog entry, as
  described in the `Usage`_ section below.
 


Implementation Details
----------------------
There are three key objects that are used to provide the ObjectCatalog functionality.

CatalogInterface
^^^^^^^^^^^^^^^^
The ``CatalogInterface`` class provides the base definitions and interface for the 
ObjectCatalog concept.
It is templated on the common base class of all derived objects that are 
creatable by the "ObjectCatalog".
In addition, ``CatalogInterface`` is templated on a variadic parameter pack that 
allows for an arbitrary constructor argument list as shown in the declaration shown below:

.. literalinclude:: ../ObjectCatalog.hpp
   :language: c++
   :start-after: //START_SPHINX_0
   :end-before: {

The ``CatalogInterface`` also defines the actual catalog type using the template arguments:

.. literalinclude:: ../ObjectCatalog.hpp
   :language: c++
   :start-after: //START_SPHINX_1
   :end-before: //STOP_SPHINX

The ``CatalogInterface::CatalogType`` is a ``std::unordered_map`` with a string "key" and a value 
type that is a pointer to the CatalogInterface that represents a specific combination of 
``BASETYPE`` and constructor arguments.

After from setting up and populating the catalog, which will be described in the "Usage" section, 
the only interface with the catalog will typically be when the ``Factory()`` method is called.
The definition of the method is given as:

.. literalinclude:: ../ObjectCatalog.hpp
   :language: c++
   :start-after: //START_SPHINX_2
   :end-before: //STOP_SPHINX

It can be seen that the static ``Factory`` method is simply a wrapper that calls the virtual 
``Allocate`` method on a the catalog which is returned by ``getCatalog()``.
The usage of the ``Factory`` method will be further discussed in the `Usage`_ section.

.. note::
  The method for organizing constructing new objects relies on a common constructor list between
  the derived type and the ``BASETYPE``. 
  This means that there is a single catalog for each combination of ``BASETYPE`` and the variadic 
  parameter pack representing the constructor arguments.
  In the future, we can investigate removing this restriction and allowing for construction of 
  a hierarchy of objects with an arbitrary constructor parameter list.

CatalogEntry
^^^^^^^^^^^^
The ``CatalogEntry`` class derives from ``CatalogInterface`` and adds the a ``TYPE`` template argument
to the arguments of the ``CatalogInterface``.

.. literalinclude:: ../ObjectCatalog.hpp
   :language: c++
   :start-after: //START_SPHINX_3
   :end-before: {

The ``TYPE`` template argument is the type of the object that you would like to be able to create 
with the "ObjectCatalog".
``TYPE`` must be derived from ``BASETYPE`` and have a constructor that matches the variadic parameter
pack specified in the template parameter list.
The main purpose of the ``CatalogEntry`` is to override the ``CatalogInterface::Allocate()`` virtual 
function s.t. when key is retrieved from the catalog, then it is possible to create a new ``TYPE``.
The ``CatalogEntry::Allocate()`` function is a simple creation of the underlying ``TYPE`` as shown by 
its definition:

.. literalinclude:: ../ObjectCatalog.hpp
   :language: c++
   :start-after: //START_SPHINX_4
   :end-before: //STOP_SPHINX

CatalogEntryConstructor
^^^^^^^^^^^^^^^^^^^^^^^
The ``CatalogEntryConstructor`` is a helper class that has a sole purpose of creating a 
new ``CatalogEntry`` and adding it to the catalog.
When a new ``CatalogEntryConstructor`` is created, a new ``CatalogEntry`` entry is created and
inserted into the catalog automatically.


Usage
-----

Creating A New Catalog
^^^^^^^^^^^^^^^^^^^^^^
When creating a new "ObjectCatalog", it typically is done within the context of a specific
``BASETYPE``.
A simple example of a class hierarchy in which we would like to use the "ObjectCatalog"
to use to generate new objects is given in the unit test located in ``testObjectCatalog.cpp``.

The base class for this example is defined as:

.. literalinclude:: ../../unitTests/dataRepositoryTests/testObjectCatalog.cpp
   :language: c++
   :start-after: //START_SPHINX_BASE
   :end-before: //STOP_SPHINX

There a couple of things to note in the definition of ``Base``:

* ``Base`` has a convenience alias to use in place of the fully templated ``CatalogInterface`` name.
* ``Base`` defines a ``getCatalog()`` function that returns a static instantiation of a
  ``CatalogInterface::CatalogType``.
  The ``CatalogInterface::getCatalog()`` function actually calls this function within the base
  class.
  This means that the base class actually owns the catalog, and the ``CatalogInterface`` is only
  operating on that ``Base::getCatalog()``, and that the definition of this function is required.

Adding A New Type To The Catalog
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once a ``Base`` class is defined with the required features, the next step is to add a new derived
type to the catalog defined in ``Base``.
There are three requirements for the new type to be registered in the catalog:

* The derived type must have a constructor with the arguments specified by the
  variadic parameter pack specified in the catalog.
* There must be a static function ``static string catalogName()`` that returns the
  name of the type that will be used to as keyname when it is registered ``Base``'s catalog.
* The new type must be registered with the catalog held in ``Base``. 
  To accomplish this, a convenience macro ``REGISTER_CATALOG_ENTRY()`` is provided.
  The arguments to this macro are the name type of Base, the type of the derived class,
  and then the variadic pack of constructor arguments.

A pair of of simple derived class that have the required methods are used in the unit test.

.. literalinclude:: ../../unitTests/dataRepositoryTests/testObjectCatalog.cpp
   :language: c++
   :start-after: //START_SPHINX_DERIVED1
   :end-before: //STOP_SPHINX

.. literalinclude:: ../../unitTests/dataRepositoryTests/testObjectCatalog.cpp
   :language: c++
   :start-after: //START_SPHINX_DERIVED2
   :end-before: //STOP_SPHINX

Allocating A New Object From The Catalog
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The test function in the unit test shows how to allocate a new object of one 
of the derived types from ``Factory`` method.
Note the call to ``Factory`` is scoped by ``Base::CatalogInterface``, which is 
an alias to the full templated instantiation of ``CatalogInterface``.
The arguments for ``Factory`` 

.. literalinclude:: ../../unitTests/dataRepositoryTests/testObjectCatalog.cpp
   :language: c++
   :start-after: //START_SPHINX_TEST
   :end-before: //STOP_SPHINX

The unit test creates two new objects of type ``Derived1`` and ``Derived2`` using the 
catalogs ``Factory`` method.
Then the test checks to see that the objects that were created are of the correct type.
This unit test has some extra output to screen to help with understanding of the 
sequence of events.
The result of running this test is::

    $ tests/testObjectCatalog 
    Calling constructor for CatalogEntryConstructor< Derived1 , Base , ... >
    Calling constructor for CatalogInterface< Base , ... >
    Calling constructor for CatalogEntry< Derived1 , Base , ... >
    Registered Base catalog component of derived type Derived1 where Derived1::catalogName() = derived1
    Calling constructor for CatalogEntryConstructor< Derived2 , Base , ... >
    Calling constructor for CatalogInterface< Base , ... >
    Calling constructor for CatalogEntry< Derived2 , Base , ... >
    Registered Base catalog component of derived type Derived2 where Derived2::catalogName() = derived2
    Running main() from gtest_main.cc
    [==========] Running 1 test from 1 test case.
    [----------] Global test environment set-up.
    [----------] 1 test from testObjectCatalog
    [ RUN      ] testObjectCatalog.testRegistration
    EXECUTING MAIN
    Creating type Derived1 from catalog of Base
    calling Base constructor with arguments (1 3.14)
    calling Derived1 constructor with arguments (1 3.14)
    Creating type Derived2 from catalog of Base
    calling Base constructor with arguments (1 3.14)
    calling Derived2 constructor with arguments (1 3.14)
    EXITING MAIN
    calling Derived2 destructor
    calling Base destructor
    calling Derived1 destructor
    calling Base destructor
    [       OK ] testObjectCatalog.testRegistration (0 ms)
    [----------] 1 test from testObjectCatalog (0 ms total)
    
    [----------] Global test environment tear-down
    [==========] 1 test from 1 test case ran. (0 ms total)
    [  PASSED  ] 1 test.
    Calling destructor for CatalogEntryConstructor< Derived2 , Base , ... >
    Calling destructor for CatalogEntryConstructor< Derived1 , Base , ... >
    Calling destructor for CatalogEntry< Derived2 , Base , ... >
    Calling destructor for CatalogInterface< Base , ... >
    Calling destructor for CatalogEntry< Derived1 , Base , ... >
    Calling destructor for CatalogInterface< Base , ... >

In the preceding output, it is clear that the static catalog in ``Base::getCatalog()``
is initialized prior the execution of main, and destroyed after the completion of main.
In practice, there have been no indicators of problems due to the use of a statically 
initialized/deinitialized catalog.

.. warning:
  The ``catalog`` is a static map, which means it is statically initialized and statically de-initialized.
  This results in a restriction in that no entry into the static map may depend on another object
  in the static map, nor on another static object in general.
  This is a well known issue with the use of static objects.
  However, this is generally not an issue in the "ObjectCatalog", as the catalog is populated with 
  ``CatalogEntry`` objects, which are not dependent on other static objects.
