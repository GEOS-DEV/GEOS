.. _dataRepository:

#####################################
GEOSX Data Repository
#####################################

The GEOSX dataRepository is intended to provide the building blocks for the code structure within GEOSX. 
The dataRepository provides a general capability to store arbitrary data and objects in a hierarchical 
structure, similar to a standard file sytem.
The components/classes of the data structure that a developer will require some knowledge of are:
* **Group** The building block for the hierachy structure.
            In the filesystem analog, the Group is similar to the file folder. 
            As such a Group may contain any number of sub-Groups, and Wrappers which are contained in a
            MappedVector.
* **Wrapper/WrapperBase** A "wrapper" class that is intended to encapuslate an object for storage in
                          a Group.
                          Wrapper is templated on the type of object it encapsulates to provide  strong 
                          type safety when retriving the wrapped objects.
                          As each Wrapper class instatiaiton will be a distinct type, Wrapper derives from
                          a non-tempalted WrapperBase class that defines common interface. The WrapperBase
                          type is the specific type that is stored in the MappedVector container within a 
                          Group.
* **MappedVector** Defines a container with an interface which is an aggregation of a std::map, and a 
                   std::vector. 
                   The types contained in the MappedVector are stored in an std::vector, and are accessible
                   through an integral index lookup.
                   In addition, there is a std::map<string,INDEX> that provides a string key lookup
                   capability to the container.

* **InputFlags** A strongly typed enum that defines the relationship with the Wrapper and the XML input.
                 For example, the enum value of OPTIONAL indicates that the data in the Wrapper is to be
                 set by the XML input if the name of the Wrapper is present in the XML file.
                 The enum value of REQUIRED indicates that the data in the Wrapper is required to be set by
                 the XML input and an error will be raised if the name is not found in the XML input.
                 Finally, the enum value of NONE indicates that Wrapper is not to be set by the XML input.

* **DefaultValue** A templated class that is used to store the default value of the contents of a Wrapper.
                   A DefaultValue specification is required if the InputFlag for a Wrapper is set to 
                   OPTIONAL.

* **RestartFlags** A strongly typed enum that specifies if the contents of a Wrapper should be written to/
                   read from restart files.

dataRepository::Group
=========================
The class `dataRepository::Group` is the object that forms the basis of the code structure hierarchy.
In the file structure analagy, the `Group` class may be thought of a file folder.
In such, the `Group` contains two 