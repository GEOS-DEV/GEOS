GEOSX data structures and XML input correspondence
##################################################################


In this document, you will learn how GEOSX classes interact with external information parsed from XML files, and how to add a new XML block that can be interpreted by GEOSX. relative permeability is used as an example.


GEOSX data structure overview
===============================================================



Groups : the base class of GEOSX
-----------------------------------------------

All GEOSX classes derive from a basic class called ``Group``. Group objects allow to structure all GEOSX objects as a tree-like structure. One could think of groups as *nodes* that can bear data and have multiple children. All groups are connected together in a tree structure to inheritance and control of the scope. For instance, an object of type "flow solver" would derive from a "*base* solver" class. Each Group object is at minimum equipped with the following member properties:

 - A pointer to its parent ``Group`` called ``m_parent`` (member classes are prefixed by ``m_``),
 - Its own data, stored for flexibility in an array of generic data "wrappers" called ``m_wrappers``,
 - One or many children objects called ``m_subGroups``.

The ``m_size`` member determines the number of data wrappers in the array ``m_wrappers``. The name of the Group is stored as a string in ``m_name``. This name can be seen as the object unique ID.

.. code-block:: cpp

    class Group
    {
      // [...]

      /// the parent of this group
      Group * m_parent = nullptr;

      /// the container for all wrappers
      wrapperMap m_wrappers;

      /// The container for all sub-groups
      subGroupMap m_subGroups;
      indexType m_size;             ///< The size/length wrappers in this group
      InputFlags m_input_flags;     ///< Input flag for this group
      string m_name;                ///< the repository name of this group. This
                                    ///< is the key in the parent group.
    };

*[Source: src/coreComponents/dataRepository/Group.hpp]*

**What is a Catalog and why do we need it?**

Some classes need external information to be instantiated. This information can be supplied in the XML input files, for instance. To connect the external (XML) and internal (C++) data structures, GEOSX uses a **Catalog** that lists all externally-accessible classes. These classes are listed using a string handle called a ``CatalogName``. The catalog follows the same tree-like structure as that of groups, and list all externally-accessible groups.

**What is a CatalogName?**

The ``CatalogName`` of an object is a *handle* (here, a string) to this object's class. Most of the time, the CatalogName and the C++ class name are identical. This helps make the code easier to debug and allows the XML/C++ correspondence to be evident. But strictly speaking, the CatalogName can be anything, as long as it refers uniquely to a specific class. The CatalogName must not be confused with the object's *name* (``m_name`` is a member of an object that stores the object's unique ID, not its class). You can have several objects of the same class and hence the same ``CatalogName``, but with different names (ie. unique ID): several fluid models, several solvers, etc.

**How can I add my new externally-accessible class to the Catalog?**

To make a class reachable externally, you need to add it to the Catalog. To do so, you must equip the class in question with a ``CatalogName()`` method that returns a string. In this document, we have referred to this returned string as the object's ``CatalogName``, but really, the method ``CatalogName()`` is what matters. The catalog of all externally-accessible objects is built by listing all ``CatalogName()`` return values. The catalog is not a flat list of strings that connect XML to C++ classes. The catalog follows the same tree-like structure as the class diagram, and thus the same rules of scope apply. For instance, the search in the catalog is always done inside the scope where the search was triggered (more on this later).


.. code-block:: cpp

  class cppNameOfMySolver : public FlowSolverBase  // <- THIS IS THE CLASS NAME
  {
  public:
    // [...]
    static string CatalogName() { return "xmlNameOfMySolver"; }  // <- THIS IS THE CATALOG (=XML) NAME
    // [...]
  }

.. highlights::

  Summary: All GEOSX objects form a tree-like structure. If an object needs to be accessible externally, it must be registered in a Catalog. The CatalogName() method returns a string handle to the object's class. The catalog has the same tree structure as the class diagram.


Registration: parsing XML input files to instantiate GEOSX objects
-------------------------------------------------------------------------


In this section, we describe the connection between **internal GEOSX objects** and **external XML tags** parsed from parameter files. We call this process *Registration*. The registration process works in three steps:

  #. The XML document is parsed. Each time a new XML tag is found, the Catalog currently in scope is inspected. The goal is to find a ``CatalogName`` string that matches the XML tag.
  #. If it is the case (the Catalog currently in scope contains a CatalogName identical to the XML tag), then the code creates a new instance of the class that the CatalogName refers to. This new object is inserted in the Group tree structure.
  #. By parsing the XML attributes of the tag, the new object properties are populated. Some checks are performed to ensure that the data supplied is conform, and that all required information is present.


*A note on Catalog scope:*

The catalog mirrors the Group tree structure. This means that for nested XML tags, the class linked to the inner tags must derive from the class linked to the outer tags. This provides consistency between the XML structure and the class diagrams.


Let's look at this process in more details.

Creating a new object and giving it a Catalog name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider registering a new physical solver deriving from a ``FlowSolverBase``. We will call this new  solver class ``cppNameOfMySolver``. This choice of name is not recommended (we want names that reflect what the solver does!), but for this particular example, we just need to know that this name is the class name inside the C++ code.

To specify parameters of this new solvers from inside an XML file, we need to be sure that the XML tag and the CatalogName of the class are the identical. Therefore, we equip the ``cppNameOfMySolver`` class with a ``CatalogName()`` method that returns the solver CatalogName (=XML name). Here, this method returns the string ``xmlNameOfMySolver``.

We have deliberately distinguished the class name from the catalog/XML name for the sake of clarity in this example. It is nevertheless a best practice to use the same name for the class and for the CatalogName.

.. code-block:: cpp

  class cppNameOfMySolver : public FlowSolverBase  // <- THIS IS THE CLASS NAME, DERIVING FROM FlowSolverBase
  {
  public:
    // [...]
    static string CatalogName() { return "xmlNameOfMySolver"; }  // <- THIS IS THE CATALOG (=XML) NAME
    // [...]
  }



Parsing XML and searching the Catalog in scope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that we have implemented a ``CatalogName()`` method returning a specific string, we can have a block in our XML input file with a tag that corresponds to the CatalogName ``xmlNameOfMySolver`` and designed to instantiate a class of type ``cppNameOfMySolver``:

.. code-block:: XML

  <?xml version="1.0" ?>
  <!--# # -->

  <Problem>
    <Solvers
      gravityVector="0.0, 0.0, -9.81">
      <xmlNameOfMySolver name="nameOfThisSolverInstance"
                               verboseLevel="1"
                               gravityFlag="1"
                               temperature="297.15"
        <SystemSolverParameters newtonTol="1.0e-6"
                                maxIterNewton="15"
                                useDirectSolver="1"/>
      </xmlNameOfMySolver>
    </Solvers>

When parsing the XML, the code will look for a CatalogName value "xmlNameOfMySolver" **inside the "Solvers" catalog**, which is in turn located **inside the "Problem" Catalog**. Again, that is because the XML structure mirrors the class diagram, and all catalogs have a local scope.


Instantiating the new solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have specified an XML block for the solver type ``xmlNameOfMySolver``. This object name was found in the local catalog, and it was dtermined to correspond to a ``cppNameOfMySolver`` object. Now, when reading the XML and encountering a ``xmlNameOfMySolver`` solver, we add a new instance of the class ``cppNameOfMySolver`` in the tree structure. It is registered as a child of its base solver class. To do this, the method ``CreateChild`` is used.


.. code-block:: cpp

    // Variable values in this example:
    // --------------------------------
    // childKey = "xmlNameOfMySolver" (string)
    // childName = "nameOfThisSolverInstance" (string)
    // SolverBase::CatalogInterface = the Catalog attached to the base Solver class
    // hasKeyName = bool method to test if the childKey string is present in the Catalog
    // RegisterGroup = method to create a new instance of the solver and add it to the group tree

    Group * PhysicsSolverManager::CreateChild( string const & childKey, string const & childName )
    {
      Group * rval = nullptr;
      if( SolverBase::CatalogInterface::hasKeyName(childKey) )
      {
        GEOS_LOG_RANK_0("Adding Solver of type " << childKey << ", named " << childName);
        rval = RegisterGroup( childName,
                              SolverBase::CatalogInterface::Factory( childKey, childName, this ) );
      }
      return rval;
    }


Note that several instances of the same type of solver can be created, as long as they each have a different name.


Filling the objects with data (wrappers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After finding and placing the new solver Group in the tree hierarchy, properties are read and stored. This is done by registering *data wrappers*. The method used to do that is called ``registerWrapper`` and is usually done in the class constructor. Note that some properties are registered at the current class level, and other properties can also be registered at a parent class level.

Here, the only data (=wrapper) that is defined at the level of our ``cppNameOfMySolver`` class is temperature, and everything else is registered at the base class (parent) level. We register a property of temperature, corresponding to the member class ``m_temperature`` of ``cppNameOfMySolver``. The registration also checks if a property is required or optional (here, it is required), and provides a brief description.

.. code-block:: cpp

  cppNameOfMySolver::cppNameOfMySolver( const string & name, Group * const parent )
    :
    FlowSolverBase( name, parent ),
  {
    this->registerWrapper( viewKeyStruct::temperatureString, &m_temperature, false )->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Temperature");
  }


This operation is done recursively if XML tags are nested.

To summarize:
-----------------------------------------------


  - Every classes in GEOSX derive from a Group in a tree-like structure. A group must have parents, can have data (in wrappers) and can have one or many children. Groups that are accessible externally are identified externally by their ``CatalogName``.
  - When parsing XML input files, GEOSX inspects each object's catalog to find classes with the same    ``CatalogName`` as the XML tag. Once it finds an XML tag in its catalog, it registers it inside the tree structure.
  - In the registration process, properties from the XML file are parsed and used to allocate member data wrappers and fully instantiate the Group class.
  - If XML tags are nested, sub-classes are allocated and processed in a nested manner.

The correspondence between XML and class hierarchy is thus respected, and the object hierarchy mirrors the XML structure.



Example: adding a new relative permeability model
================================================================================

This example is taken from the class ``BrooksCoreyRelativePermeability``, derived from RelativePermeabilityBase.

*[source: src/coreComponents/constitutive/RelPerm/BrooksCoreyRelativePermeability.hpp/.cpp]*

Implement a ``CatalogName`` function (.hpp):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::cpp

  static std::string CatalogName() { return dataRepository::keys::brooksCoreyRelativePermeability; }

Here the namespace is defined at the top of the hpp file as:

.. code-block::cpp

  namespace dataRepository
  {
  namespace keys
  {
  string const brooksCoreyRelativePermeability = "BrooksCoreyRelativePermeability";
  }
  }

Now every time a "BrooksCoreyRelativePermeability" string is encountered inside a Relative Permeability catalog, we will instantiate a class ``BrooksCoreyRelativePermeability``.


Register all data wrappers (.hpp):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When attaching properties (ie. data wrappers) to a class, a similar registration process must be done. Every property is accessed through its ``ViewKey`` name space.

.. code-block::cpp

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString = "phaseMinVolumeFraction";  // = XML attribute
    static constexpr auto phaseRelPermExponentString   = "phaseRelPermExponent";  // = XML attribute
    static constexpr auto phaseRelPermMaxValueString   = "phaseRelPermMaxValue";  // = XML attribute

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseMinVolumeFraction = { phaseMinVolumeFractionString };  // = struct key name
    ViewKey phaseRelPermExponent   = { phaseRelPermExponentString };  // = struct key name
    ViewKey phaseRelPermMaxValue   = { phaseRelPermMaxValueString };  // = struct key name

  } vieKeysBrooksCoreyRelativePermeability;


Define data members (.hpp):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The data members are defined in the class.

.. code-block::cpp

  protected:
    array1d<real64> m_phaseMinVolumeFraction;
    array1d<real64> m_phaseRelPermMaxValue;
    array1d<real64> m_phaseRelPermExponent;
    real64 m_satScale;


Implement the data registration process (registerWrapper):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The registration process connects the attributes values in the XML file to class member data:


.. code-block:: cpp

  BrooksCoreyRelativePermeability::BrooksCoreyRelativePermeability( std::string const & name, Group * const parent )
    : RelativePermeabilityBase( name, parent )
  {
    registerWrapper( viewKeyStruct::phaseMinVolumeFractionString, &m_phaseMinVolumeFraction, false )->
      setApplyDefaultValue(0.0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Minimum volume fraction value for each phase");

    registerWrapper( viewKeyStruct::phaseRelPermExponentString,   &m_phaseRelPermExponent,   false )->
      setApplyDefaultValue(1.0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("MinimumRel perm power law exponent for each phase");

    registerWrapper( viewKeyStruct::phaseRelPermMaxValueString,   &m_phaseRelPermMaxValue,   false )->
      setApplyDefaultValue(0.0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Maximum rel perm value for each phase");
  }


*For more examples on how to contribute to GEOSX, please read "Adding a new Physics Solver".*
