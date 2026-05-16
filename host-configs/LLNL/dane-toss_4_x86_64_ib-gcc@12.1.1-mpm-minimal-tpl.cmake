# Dane toss_4_x86_64_ib gcc@12.1.1 MPM Minimal-TPL host-config.
# This starts from the generated Dane host-config and then strips non-MPM packages/TPLs.

include(${CMAKE_CURRENT_LIST_DIR}/dane-toss_4_x86_64_ib-gcc@12.1.1.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../profiles/mpm-minimal-tpl.cmake)
