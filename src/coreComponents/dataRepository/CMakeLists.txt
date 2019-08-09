
#
# Specify all headers
#
set(dataRepository_headers
    BufferOps_inline.hpp
    BufferOps.hpp
    ExecutableGroup.hpp
    Group.hpp
    InputFlags.hpp
    KeyIndexT.hpp
    MappedVector.hpp
    ReferenceWrapper.hpp
    RestartFlags.hpp
    SidreWrapper.hpp
    Wrapper.hpp
    WrapperBase.hpp
    xmlWrapper.hpp
    )

#
# Specify all sources
#
set(dataRepository_sources
    Group.cpp
    SidreWrapper.cpp
    Wrapper.cpp
    WrapperBase.cpp
    xmlWrapper.cpp
    )

set( dependencyList common codingUtilities pugixml )

blt_add_library( NAME                  dataRepository
                 SOURCES               ${dataRepository_sources}
                 HEADERS               ${dataRepository_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                TRUE
               )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX dataRepository )

add_subdirectory( unitTests )
