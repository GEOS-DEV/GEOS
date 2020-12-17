# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Adiak(CMakePackage):
    """Adiak collects metadata about HPC application runs and provides it
       to tools."""

    homepage = "https://github.com/LLNL/Adiak"
    git      = "https://github.com/LLNL/Adiak.git"

    version('master', branch='master', submodules='True')
    version('0.2.1', tag='v0.2.1', submodules="True")
    version('0.2.0', tag='v0.2.0', submodules="True")
    version('0.1.1', tag='v0.1', submodules="True")

    variant('mpi', default=True, description='Build with MPI support')
    variant('shared', default=True, description='Build dynamic libraries')

    depends_on('mpi', when='+mpi')

    def cmake_args(self):
        args = []
        if self.spec.satisfies('+mpi'):
            args.append('-DMPICXX=%s' % self.spec['mpi'].mpicxx)
            args.append('-DMPICC=%s' % self.spec['mpi'].mpicc)
            args.append('-DENABLE_MPI=ON')
        else:
            args.append('-DENABLE_MPI=OFF')

        if (self.spec.satisfies('+shared')):
            args.append('-DBUILD_SHARED_LIBS=ON')
        else:
            args.append('-DBUILD_SHARED_LIBS=OFF')

        args.append('-DENABLE_TESTS=OFF')
        return args
