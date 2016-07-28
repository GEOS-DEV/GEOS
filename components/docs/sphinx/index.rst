
Sidre User Documentation
=========================

The Sidre (Simulation data repository) component of the CS Toolkit provides
tools to centralize data management in HPC applications (i.e., description, 
allocation, access, etc.). Requirements for Sidre grew out of related
functionality that is central to several current LLNL applications as well 
needs identified for new codes targeted for future advanced architectures 
that must manage data allocation and placement carefully to run efficiently on
those platforms. Unlike capabilities in existing codes, which were typically 
developed specifically for each code, Sidre is designed from inception to be 
shared by different applications. The goal of Sidre is to facilitate 
coordination and sharing of data: across physics packages in integrated 
applications, and between applications, libraries, and tools that provide 
capabilities such as file I/O, in situ visualization and analysis, etc.


Introduction
-------------

Sidre provides simple application-level semantics to describe, 
allocate/deallocate, and provide access to data. Currently supported
capbilities include:

* Separate data description and allocation operations. This allows applications
  to describe data they need and then perform "counting exercises" to decide
  how best to place the data in memory. 
* Externally-owned "opaque" or described data. Sidre can accept a pointer to 
  an externally-allocated data object and provide access to it by name. When
  the external data is described to Sidre, interaction with it by the user
  is similar in many ways to data that Sidre owns. When the data is not 
  described (i.e., it is "opaque"), Sidre can provide access to the data 
  via a pointer, but the consumer of the pointer must know type information
  to do anything substantial with the data.
* Multiple, different "views" into a chunk of (shared) data. A Sidre view 
  includes description semantics to define data type, number of elements, 
  offset, stride, etc. Thus, a chunk of data in memory can be interpreted 
  conceptually in different ways.
* Arbitrarily-organized data hierarchies. Many mesh-based application codes 
  organize data into hierarchies of contexts (e.g., domains, regions, blocks, 
  mesh centerings, subsets of elements containing different 
  materials, etc.). Sidre supports hierarchical, tree-based organizations 
  in a simple, flexible way.
* APIs for C++, C, and Fortran along with mechanisms to insure inter-language
  data consistency.

So far, Sidre development has focused on designing and building flexible and
powerful concepts to build on. The Sidre API includes four main concepts:

* **Datastore.**  This is main interface. It is small and simple.
* **Buffer.**  Describes and holds a chunk of data in memory owned by Sidre.
* **View.**   Describes a conceptual layout of data in memory.
* **Group.**  Used to define a hierachical collection of data. A Group can have one parent Group and one or more child Groups. A Group can also have a collection of Views.

These concepts will be described in more detail in later sections.

At this point, Sidre supports simple data types such as scalars, strings, and 
(multidimensional) arrays of scalars. Fundamentally, Sidre does not preclude 
the use of more complex data structures, but does not directly support them
currently. Future Sidre development will expand core functionality to support 
additional features such as:

* Mechanisms to associate data with memory spaces and tranfer data between
  spaces.
* Support for more complex data types.
* "Attributes" whereby data objects can be tagged with a key (e.g., string) 
  and an optional value. This would provide a potentially more convenient 
  mechanism for applications to retrieve and perform operations on a 
  collection of data objects.

Support for these enhancements and others will be added based on application
needs and use cases.


**Contents:**

.. toctree::
   :maxdepth: 2

   first_example
   core_concepts
   file_io
   sidre_conduit

