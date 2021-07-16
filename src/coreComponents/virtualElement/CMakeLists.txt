#
# Specify all headers
#
set( virtualElement_headers
     ConformingVirtualElementOrder1.hpp
     ConformingVirtualElementOrder1_impl.hpp
     VirtualElementBase.hpp
   )

set( dependencyList mesh )

blt_add_library( NAME        virtualElement
                 HEADERS     ${virtualElement_headers}
                 DEPENDS_ON  ${dependencyList}
               )

geosx_add_code_checks( PREFIX virtualElement )
