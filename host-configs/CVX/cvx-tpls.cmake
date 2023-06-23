# Pavel: edit that to provide the path for thirdPartyLibs
#        alternatively (and preferably) use
#          python scripts/config-build.py -hc host-configs/cvx-wsl.cmake -bt Release -D GEOSX_TPL_DIR=/full/path/to/thirdPartyLibs
#set(GEOSX_TPL_DIR "/home/ptls/thirdPartyLibs/install-cvx-wsl-debug" CACHE PATH "" FORCE)

include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
