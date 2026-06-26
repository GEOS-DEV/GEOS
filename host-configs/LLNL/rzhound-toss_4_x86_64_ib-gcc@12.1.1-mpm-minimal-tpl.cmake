# RZHound toss_4_x86_64_ib gcc@12.1.1 MPM Minimal-TPL host-config.
#
# RZHound is a CPU-only TOSS4 Sapphire Rapids Slurm cluster.  For the current
# MPM minimal-TPL CPU build, reuse the Dane CPU host-config wrapper.  Keep this
# file separate so setupMPM can choose a RZHound-named host-config and build
# directory while still sharing the same GEOS options.

include(${CMAKE_CURRENT_LIST_DIR}/dane-toss_4_x86_64_ib-gcc@12.1.1-mpm-minimal-tpl.cmake)
