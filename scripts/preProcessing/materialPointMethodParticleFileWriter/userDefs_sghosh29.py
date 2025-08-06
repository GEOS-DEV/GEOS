# -*- coding: utf-8 -*-

# Don't modify this file.
# Make a copy of this file renamed with your LC username.  Don't add your version to the git repo.
import platform
lassen = 'lassen' in platform.node()

#Should the particleStressHistory be in a sbatch wrapper so if the os hangs it isn't using job time
  
if lassen:
  geosPath='/home/sghosh29/data-rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx'
  testRunDirectory='/home/sghosh29/data-rhurley6/sohanjit/testGEOS/'
  defaultRunDirectory='/home/sghosh29/data-rhurley6/sohanjit/testGEOS/'
  defaultBank='cbronze'
else:
  geosPath='/home/sghosh29/data-rhurley6/sohanjit/GEOS/build-rockfish-clang-release/bin/geosx'
  testRunDirectory='/home/sghosh29/data-rhurley6/sohanjit/testGEOS/'
  defaultRunDirectory='/home/sghosh29/data-rhurley6/sohanjit/testGEOS/'
  defaultBank='imcomp'
