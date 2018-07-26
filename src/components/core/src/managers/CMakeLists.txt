
#
# Specify all headers
#
set( managers_headers
     DomainPartition.hpp
     EventManager.hpp
     Functions/FunctionBase.hpp
     Functions/SymbolicFunction.hpp
     Functions/TableFunction.hpp
     Functions/CompositeFunction.hpp
	 Functions/NewFunctionManager.hpp
     ObjectManagerBase.hpp
     ProblemManager.hpp
	 NumericalMethodsManager.hpp
   )


#
# Specify all sources
#
set( managers_sources
     DomainPartition.cpp
     EventManager.cpp
     Functions/FunctionBase.cpp
     Functions/SymbolicFunction.cpp
     Functions/TableFunction.cpp
     Functions/CompositeFunction.cpp
     Functions/NewFunctionManager.cpp
     ObjectManagerBase.cpp
     ProblemManager.cpp
	 NumericalMethodsManager.cpp
   )

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
               
target_include_directories( managers PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src )

geosx_add_code_checks( PREFIX managers )