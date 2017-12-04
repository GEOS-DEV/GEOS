#!/bin/bash
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-clang\@4.0.0.cmake -tpl -bt Release -ip /usr/gapps/GEOS/geosx/thirdPartyLibs/install-toss_3_x86_64_ib-clang\@4.0.0-release
python scripts/config-build.py -hc host-configs/toss_3_x86_64_ib-gcc\@7.1.0.cmake   -tpl -bt Release -ip /usr/gapps/GEOS/geosx/thirdPartyLibs/install-toss_3_x86_64_ib-gcc\@7.1.0-release

