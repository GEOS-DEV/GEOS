# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Essl(Package):
    """IBM's Engineering and Scientific Subroutine Library (ESSL)."""

    # homepage = "https://www.ibm.com/systems/power/software/essl/"

    variant('ilp64', default=False, description='64 bit integers')
    variant(
        'threads', default='openmp',
        description='Multithreading support',
        values=('openmp', 'none'),
        multi=False
    )
    variant('cuda', default=False, description='CUDA acceleration')
    variant('lapack', default=False, description='Support lapack interface.')

    depends_on('cuda', when='+cuda')

    provides('blas')
    provides('lapack', when='+lapack')

    conflicts('+cuda', when='+ilp64',
              msg='ESSL+cuda+ilp64 cannot combine CUDA acceleration'
              ' 64 bit integers')

    conflicts('+cuda', when='threads=none',
              msg='ESSL+cuda threads=none cannot combine CUDA acceleration'
              ' without multithreading support')

    @property
    def blas_libs(self):
        spec = self.spec

        if '+ilp64' in spec:
            essl_lib = ['libessl6464']
        else:
            essl_lib = ['libessl']

        if spec.satisfies('threads=openmp'):
            # ESSL SMP support requires XL or Clang OpenMP library
            if '%xl' in spec or '%xl_r' in spec or '%clang' in spec:
                if '+ilp64' in spec:
                    essl_lib = ['libesslsmp6464']
                else:
                    if '+cuda' in spec:
                        essl_lib = ['libesslsmpcuda']
                    else:
                        essl_lib = ['libesslsmp']

        essl_libs = find_libraries(
            essl_lib,
            root=self.prefix.lib64,
            shared=True
        )

        if '+cuda' in spec:
            essl_libs += spec['cuda'].libs

        raise RuntimeError("%s" % essl_libs)
        return essl_libs

    @property
    def lapack_libs(self):
        essl_libs = self.blas_libs
        essl_libs += find_libraries(
              ['liblapackforessl'],
              root=self.prefix.lib64,
              shared=True
        )

        return essl_libs
    
    @property
    def libs(self):
        result = self.blas_libs
        if '+lapack' in self.spec:
            result += self.lapack_libs

        return result

    def install(self, spec, prefix):
        raise InstallError('IBM ESSL is not installable;'
                           ' it is vendor supplied')
