
#
# Specify all headers
#
set(managers_headers
    DomainPartition.hpp
    EventManager.hpp
    Events/EventBase.hpp
    Events/PeriodicEvent.hpp
    Events/HaltEvent.hpp
    Events/SoloEvent.hpp
    Outputs/OutputManager.hpp
    Outputs/OutputBase.hpp
    Outputs/SiloOutput.hpp
    Outputs/RestartOutput.hpp
    Outputs/BlueprintOutput.hpp
    Outputs/ChomboIO.hpp
    Functions/FunctionBase.hpp
    Functions/SymbolicFunction.hpp
    Functions/TableFunction.hpp
    Functions/CompositeFunction.hpp
    Functions/NewFunctionManager.hpp
    ObjectManagerBase.hpp
    ProblemManager.hpp
    NumericalMethodsManager.hpp
    FieldSpecification/FieldSpecificationBase.hpp
    FieldSpecification/FieldSpecificationManager.hpp
    FieldSpecification/DirichletBoundaryCondition.hpp)


#
# Specify all sources
#
set(managers_sources
  	DomainPartition.cpp
    EventManager.cpp
    Events/EventBase.cpp
    Events/PeriodicEvent.cpp
    Events/HaltEvent.cpp
    Events/SoloEvent.cpp
    Outputs/OutputManager.cpp
    Outputs/OutputBase.cpp
    Outputs/SiloOutput.cpp
    Outputs/RestartOutput.cpp
    Outputs/BlueprintOutput.cpp
    Outputs/ChomboIO.cpp
    Functions/FunctionBase.cpp
    Functions/SymbolicFunction.cpp
    Functions/TableFunction.cpp
    Functions/CompositeFunction.cpp
    Functions/NewFunctionManager.cpp
    ObjectManagerBase.cpp
    ProblemManager.cpp
    NumericalMethodsManager.cpp
    FieldSpecification/FieldSpecificationBase.cpp
    FieldSpecification/FieldSpecificationManager.cpp
    FieldSpecification/DirichletBoundaryCondition.cpp)

#set( dependencyList dataRepository fileIO  )
#set( dependencyList2 optionparser mathpresso RAJA trilinos )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )


blt_add_library( NAME                  managers
                 SOURCES               ${managers_sources}
                 HEADERS               ${managers_headers}
                 DEPENDS_ON            dataRepository fileIO optionparser mathpresso RAJA trilinos
#                 SHARED FALSE
               )


#get_target_property( id managers INCLUDE_DIRECTORIES )
#message("managers INCLUDE_DIRECTORIES = ${id}")
#
#get_target_property( iid managers INTERFACE_INCLUDE_DIRECTORIES )
#message("managers INTERFACE_INCLUDE_DIRECTORIES = ${iid}")
#
#get_target_property( isid managers INTERFACE_SYSTEM_INCLUDE_DIRECTORIES )
#message("managers INTERFACE_SYSTEM_INCLUDE_DIRECTORIES = ${isid}")
#
#
#
#get_target_property( ill managers INTERFACE_LINK_LIBRARIES )
#message("managers INTERFACE_LINK_LIBRARIES = ${ill}")
#
#foreach( dependency ${ill} )
#    get_target_property( depList ${dependency} INTERFACE_INCLUDE_DIRECTORIES )
#    message( "    ${dependency} iid = ${depList} ")
#    list( APPEND _expanded_DEPENDS_ON ${depList} )
#endforeach()
 



target_include_directories( managers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


geosx_add_code_checks( PREFIX managers )
