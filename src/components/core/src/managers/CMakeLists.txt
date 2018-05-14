
#
# Specify all headers
#
set(managers_headers
    DomainPartition.hpp
    EventManager.hpp
    Functions/FunctionBase.hpp
    Functions/SymbolicFunction.hpp
    Functions/TableFunction.hpp
    Functions/CompositeFunction.hpp
    Functions/NewFunctionManager.hpp
    ObjectManagerBase.hpp
    ProblemManager.hpp )


#
# Specify all sources
#
set(managers_sources
  	DomainPartition.cpp
    EventManager.cpp
    Functions/FunctionBase.cpp
    Functions/SymbolicFunction.cpp
    Functions/TableFunction.cpp
    Functions/CompositeFunction.cpp
    Functions/NewFunctionManager.cpp
    ObjectManagerBase.cpp
    ProblemManager.cpp )




blt_add_library( NAME                  managers
                 SOURCES               ${managers_sources}
                 HEADERS               ${managers_headers}
                 DEPENDS_ON            dataRepository fileIO optionparser mathpresso RAJA trilinos
                 SHARED FALSE
               )