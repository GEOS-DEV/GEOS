#!/bin/bash
python scripts/config-build.py -hc host-configs/chaos_5_x86_64_ib-gcc@4.9.3.cmake -tpl
python scripts/config-build.py -hc host-configs/chaos_5_x86_64_ib-gcc@4.9.3.cmake

python scripts/config-build.py -hc host-configs/chaos_5_x86_64_ib-gcc@4.9.3.cmake -tpl -bt Release
python scripts/config-build.py -hc host-configs/chaos_5_x86_64_ib-gcc@4.9.3.cmake -bt Release