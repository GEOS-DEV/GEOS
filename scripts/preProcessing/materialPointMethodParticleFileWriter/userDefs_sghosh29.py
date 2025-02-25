# -*- coding: utf-8 -*-

# Don't modify this file.
# Make a copy of this file renamed with your LC username.  Don't add your version to the git repo.
import platform
lassen = 'lassen' in platform.node()

#Should the particleStressHistory be in a sbatch wrapper so if the os hangs it isn't using job time
  
if lassen:
  geosPath='/data1/sghosh29/Working_MPM_LLNL/GEOS/build-system76-pc-clang-release/bin/geosx'
  testRunDirectory='/data1/sghosh29/Working_MPM_LLNL/testGEOS'
  defaultRunDirectory='/data1/sghosh29/Working_MPM_LLNL/testGEOS'
  defaultBank='cbronze'
else:
  geosPath='/data1/sghosh29/Working_MPM_LLNL/GEOS/build-system76-pc-clang-release/bin/geosx'
  testRunDirectory='/data1/sghosh29/Working_MPM_LLNL/testGEOS'
  defaultRunDirectory='/data1/sghosh29/Working_MPM_LLNL/testGEOS'
  defaultBank='imcomp'
