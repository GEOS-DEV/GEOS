#
# Specify all headers
#
set( finiteElement_sources
     FiniteElementSpace.hpp
     FiniteElementManager.hpp
     ElementLibrary/FiniteElement.h
     ElementLibrary/FiniteElementBase.h
     Kinematics.h
     basis/BasisBase.hpp
     basis/LagrangeBasis.hpp
     basis/Polynomial.hpp    
     quadrature/QuadratureBase.hpp
     quadrature/GaussQuadrature.hpp
     FiniteElementSpaceManager.hpp
     basis/BasisFunctionManager.hpp
     quadrature/QuadratureRuleManager.hpp
   )
#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementSpace.cpp
     FiniteElementManager.cpp
     ElementLibrary/FiniteElement.cpp
     ElementLibrary/FiniteElementBase.cpp
     Kinematics.cpp
     basis/BasisBase.cpp
     basis/LagrangeBasis.cpp
     basis/Polynomial.cpp    
     quadrature/QuadratureBase.cpp
     quadrature/GaussQuadrature.cpp
     FiniteElementSpaceManager.cpp
     basis/BasisFunctionManager.cpp
     quadrature/QuadratureRuleManager.cpp
   )

#set( dependencyList dataRepository codingUtilities )
#set( dependencyList2 "" )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )

blt_add_library( NAME                  finiteElement
                 SOURCES               ${finiteElement_sources}
                 HEADERS               ${finiteElement_headers}
                 DEPENDS_ON            dataRepository codingUtilities
                 SHARED FALSE
               )
               
target_include_directories( finiteElement PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX finiteElement )