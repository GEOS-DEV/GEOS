###############################################################################
GEOSX Developer Guildines
###############################################################################



Language
==================
GEOSX is written in standard c++14.
We will consider c++17 features, depending on compiler support on LLNL systems.

Target Platforms
================
- Linux
- BlueOS
- MacOSX
- Windows???

Naming Conventions
=================================

File Names
----------------
- C++ header files are always named with a file extion of  \*.hpp.
- C++ header implementation files, which contain templated or inline function definitions, are always named \*_impl.hpp.
- C++ source files are always named with a file extion of  \*.cpp.
- C++ class declarations and defintions are contained files with identical names, except for the extensions.
- C++ free function headers and source files are declared/defined in files with identical names, except for the extension.

For example, a class named "Foo" may be declared in a file named "Foo.hpp", with inline/templated functions 
defined in "Foo_impl.hpp, with the source implmenetaion contained in Foo.cpp.

- There should not be identical filenames that only differ by case. Some filesystems are not case-sensitive, 
  and worse, some filesystems such as MacOSX are case-preserving but not case sensitive.

Function Names
--------------
Function and member function names should be `PascalCase <https://en.wikipedia.org/wiki/Camel_case/>`__ 
or `snake_case <https://en.wikipedia.org/wiki/Snake_case/>`__.
Please be consistent with your choice within your files.

Variable Names
--------------
Variables should be either `camelCase <https://en.wikipedia.org/wiki/Camel_case/>`__, or 
`snake_case <https://en.wikipedia.org/wiki/Snake_case/>`__."
Whatever you choose, please be consistent with your choice in the functional scope you are working in.

Member data should be `camelCase <https://en.wikipedia.org/wiki/Camel_case/>`__ prefix with 'm\_' (i.e. double m_dataVariable;)

Type Names
--------------
Please use `PascalCase <https://en.wikipedia.org/wiki/Camel_case/>`__ for typenames (i.e. classes)

   .. code-block:: c
   
      class MyClass;



   .. code-block:: c
   
      class MyClass
      {
        double m_doubleDataMember;
        int m_integerDataMember;
      }

Namespace Names
----------------
Namespaces names are all lower `camel case <https://en.wikipedia.org/wiki/Camel_case/>`__.

Example
-------
One example of would be a for a class named "Foo", the declaration would be in a header file named "Foo.hpp"

.. code-block:: c

  /*
   * Foo.hpp
   */

  namespace bar
  {
  
  class Foo
  {
  public:
    Foo();
  private:
    double m_myDouble;
  }
  }

and a source file named "Foo.cpp"

.. code-block:: c

  /*
   * Foo.cpp
   */
  namespace bar
  {
    Foo::Foo():
      m_myDouble(0.0)
    {
      // some constructor stuff
    }
  }

Code Format
=================================
GEOSX applies a variant of the 
`BSD/Allman Style <https://en.wikipedia.org/wiki/Indentation_style#Allman_style/>`__.
Key points to the GEOSX style are:

#. Opening braces ( i.e. { ) go on the next line of any control statment, and are not indented from the control statement .
#. No Tabs. Only spaces.
#. 2-space indentation

   .. code-block:: c

      for( int i=0 ; i<10 ; ++i )
      {
        std::cout<<"blah"<<std::endl;
      } 

#. Try to stay under 100 character line lengths. To achieve this apply these rules in order
#. Align function declaration/definitions/calls on argument list
#. Break up return type and function definition on new line
#. Break up scope resolution operators

   .. code-block:: c
    
    void 
    SolidMechanics_LagrangianFEM::
    TimeStepExplicit( real64 const& time_n,
                      real64 const& dt,
                      const int cycleNumber,
                      DomainPartition * const domain )
     {
       code here 
     }

As part of the continuous integration testing, this GEOSX code style is enforced via the uncrustify tool.
While quite extensive, uncrustify does not enforce every example of the preferred code style. 
In cases where uncrusitfy is unable to enforce code style, it will ignore formatting rules.
In these cases it is acceptible to proceed with pull requests, as there is no logical recourse.