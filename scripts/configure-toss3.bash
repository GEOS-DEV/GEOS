#!/bin/bash
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-clang\@4.0.0.cmake
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-clang\@4.0.0.cmake -bt Release
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-clang\@4.0.0-NoOPENMP.cmake
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-clang\@4.0.0-NoOPENMP.cmake -bt Release
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-clang\@6.0.0.cmake
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-clang\@6.0.0.cmake -bt Release
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-gcc\@8.1.0.cmake
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-gcc\@8.1.0.cmake   -bt Release

