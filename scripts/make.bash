#!/bin/bash
python scripts/config-build.py -hc host-configs/darwin-clang.cmake
python scripts/config-build.py -hc host-configs/darwin-clang37.cmake
python scripts/config-build.py -hc host-configs/darwin-gcc5.cmake
python scripts/config-build.py -bp build-xcode -hc host-configs/darwin-clang.cmake -x
