.. _dataRepository:

#####################################
GEOSX Data Repository
#####################################

The GEOSX "dataRepository" is intended to provide the building blocks for the code structure within GEOSX. 
The "dataRepository" provides a general capability to store arbitrary data and objects in a hierarchical 
structure, similar to a standard file system.
The components/classes of the data structure that a developer will require some knowledge of are:

* :ref:`Group`
  The inheritable building block for the hierarchy. 
  In the filesystem analog, the Group is similar to the file folder.
  
* **Wrapper/WrapperBase** 
  A class that is intended to encapsulate an object for storage in a Group, as well
  as providing an interface for performing common operations on that object.
  Wrapper is templated on the type of object it encapsulates, thus providing strong 
  type safety when retrieving the objects.
  As each Wrapper class instantiation will be a distinct type, Wrapper derives from
  a non-templated WrapperBase class that defines a common interface. 
  The WrapperBase type is the type of pointer that is stored in the MappedVector 
  container within a Group.
                          
* **MappedVector** 
  Defines a container that provides an interface which is an aggregation of a
  ``std::unordered_map`` and a ``std::vector``.
  The types contained in the MappedVector are stored in a ``std::vector``, and thus are 
  accessible through an integral index lookup.
  In addition, there is a ``std::unordered_map<string,INDEX>`` member that provides a string 
  key lookup capability to the container if that is the preferred interface.
  Note that the lookup will have some overhead assocaited with the unordered_map hash 
  lookup.

* **InputFlags** 
  A strongly typed enum that defines the relationship with the Wrapper and the XML input.
  For example, the enum value of OPTIONAL indicates that the data in the Wrapper is to be
  set by the XML input if the name of the Wrapper is present in the XML file.
  The enum value of REQUIRED indicates that the data in the Wrapper is required to be set by
  the XML input and an error will be raised if the name is not found in the XML input.
  Finally, the enum value of NONE indicates that Wrapper is not to be set by the XML input.

* **DefaultValue** 
  A templated type that is used to store the default value of a Wrapper.
  The application of DefaultValue is 
  The specification of a DefaultValue is required if the InputFlag for a Wrapper is set to 
  OPTIONAL.

* **RestartFlags** 
  A strongly typed enum that specifies if the contents of a Wrapper should be written 
  to and/or read from restart files.
                   
* :ref:`ObjectCatalog` 
  A statically initialized factory used for creating new objects in the dataRepository.

.. toctree::
   :maxdepth: 1

   Group
   ObjectCatalog
