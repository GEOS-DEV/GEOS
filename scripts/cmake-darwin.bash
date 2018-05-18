#!/bin/bash
#python scripts/config-build.py -hc host-configs/darwin-clang.cmake  -tpl -bt Release
#python scripts/config-build.py -hc host-configs/darwin-clang5.cmake -tpl -bt Release
#python scripts/config-build.py -hc host-configs/darwin-gcc7.cmake   -tpl -bt Release



python scripts/config-build.py -hc host-configs/darwin-clang.cmake
python scripts/config-build.py -hc host-configs/darwin-clang5.cmake
python scripts/config-build.py -hc host-configs/darwin-gcc7.cmake
python scripts/config-build.py -bp build-xcode -hc host-configs/darwin-clang.cmake -x
