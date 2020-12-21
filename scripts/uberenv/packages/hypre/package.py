# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import os
import sys


class Hypre(Package):
    """Hypre is a library of high performance preconditioners that
       features parallel multigrid methods for both structured and
       unstructured grid problems."""

    homepage = "http://computing.llnl.gov/project/linear_solvers/software.php"
    url      = "https://github.com/hypre-space/hypre/archive/v2.14.0.tar.gz"
    git      = "https://github.com/hypre-space/hypre.git"

    maintainers = ['ulrikeyang', 'osborn9', 'balay']

    version('develop', branch='master')
    version('2.20.0', sha256='5be77b28ddf945c92cde4b52a272d16fb5e9a7dc05e714fc5765948cba802c01')
    version('2.19.0', sha256='466b19d8a86c69989a237f6f03f20d35c0c63a818776d2cd071b0a084cffeba5')
    version('2.18.2', sha256='28007b5b584eaf9397f933032d8367788707a2d356d78e47b99e551ab10cc76a')
    version('2.18.1', sha256='220f9c4ad024e815add8dad8950eaa2d8f4f231104788cf2a3c5d9da8f94ba6e')
    version('2.18.0', sha256='62591ac69f9cc9728bd6d952b65bcadd2dfe52b521081612609804a413f49b07')
    version('2.17.0', sha256='4674f938743aa29eb4d775211b13b089b9de84bfe5e9ea00c7d8782ed84a46d7')
    version('2.16.0', sha256='33f8a27041e697343b820d0426e74694670f955e21bbf3fcb07ee95b22c59e90')
    version('2.15.1', sha256='50d0c0c86b4baad227aa9bdfda4297acafc64c3c7256c27351f8bae1ae6f2402')
    version('2.15.0', sha256='2d597472b473964210ca9368b2cb027510fff4fa2193a8c04445e2ed4ff63045')
    version('2.14.0', sha256='705a0c67c68936bb011c50e7aa8d7d8b9693665a9709b584275ec3782e03ee8c')
    version('2.13.0', sha256='3979602689c3b6e491c7cf4b219cfe96df5a6cd69a5302ccaa8a95ab19064bad')
    version('2.12.1', sha256='824841a60b14167a0051b68fdb4e362e0207282348128c9d0ca0fd2c9848785c')
    version('2.11.2', sha256='25b6c1226411593f71bb5cf3891431afaa8c3fd487bdfe4faeeb55c6fdfb269e')
    version('2.11.1', sha256='6bb2ff565ff694596d0e94d0a75f0c3a2cd6715b8b7652bc71feb8698554db93')
    version('2.10.1', sha256='a4a9df645ebdc11e86221b794b276d1e17974887ead161d5050aaf0b43bb183a')
    version('2.10.0b', sha256='b55dbdc692afe5a00490d1ea1c38dd908dae244f7bdd7faaf711680059824c11')
    version('xsdk-0.2.0', tag='xsdk-0.2.0', git='https://github.com/LLNL/hypre.git')

    # Versions 2.13.0 and later can be patched to build shared
    # libraries on Darwin; the patch for this capability does not
    # apply to version 2.12.1 and earlier due to changes in the build system
    # between versions 2.12.1 and 2.13.0.
    variant('shared', default=(sys.platform != 'darwin'),
            description="Build shared library (disables static library)")
    # Use internal SuperLU routines for FEI - version 2.12.1 and below
    variant('internal-superlu', default=False,
            description="Use internal SuperLU routines")
    variant('superlu-dist', default=False,
            description='Activates support for SuperLU_Dist library')
    variant('int64', default=False,
            description="Use 64bit integers")
    variant('mixedint', default=False,
            description="Use 64bit integers while reducing memory use")
    variant('complex', default=False, description='Use complex values')
    variant('mpi', default=True, description='Enable MPI support')
    variant('openmp', default=False, description='Enable OpenMP support')
    variant('debug', default=False,
            description='Build debug instead of optimized version')

    # Patch to add ppc64le in config.guess
    patch('ibm-ppc64le.patch', when='@:2.11.1')

    # Patch to build shared libraries on Darwin
    patch('darwin-shared-libs-for-hypre-2.13.0.patch', when='+shared@2.13.0 platform=darwin')
    patch('darwin-shared-libs-for-hypre-2.14.0.patch', when='+shared@2.14.0 platform=darwin')
    patch('superlu-dist-link-2.15.0.patch', when='+superlu-dist @2.15:2.16.0')
    patch('superlu-dist-link-2.14.0.patch', when='+superlu-dist @:2.14.0')
    patch('hypre21800-compat.patch', when='@2.18.0')

    depends_on("mpi", when='+mpi')
    depends_on("blas")
    depends_on("lapack")
    depends_on('superlu-dist', when='+superlu-dist+mpi')

    # Patch to build shared libraries on Darwin does not apply to
    # versions before 2.13.0
    conflicts("+shared@:2.12.99 platform=darwin")

    # Conflicts
    # Option added in v2.13.0
    conflicts('+superlu-dist', when='@:2.12.99')

    # Internal SuperLU Option removed in v2.13.0
    conflicts('+internal-superlu', when='@2.13.0:')

    # Option added in v2.16.0
    conflicts('+mixedint', when='@:2.15.99')

    def url_for_version(self, version):
        if version >= Version('2.12.0'):
            url = 'https://github.com/hypre-space/hypre/archive/v{0}.tar.gz'
        else:
            url = 'http://computing.llnl.gov/project/linear_solvers/download/hypre-{0}.tar.gz'

        return url.format(version)

    def install(self, spec, prefix):
        # Note: --with-(lapack|blas)_libs= needs space separated list of names
        lapack_libs = spec['lapack'].libs
        blas_libs = spec['blas'].libs

        configure_args = [
            '--prefix=%s' % prefix,
            '--with-lapack-libs=%s' % ' '.join(lapack_libs.names),
            '--with-lapack-lib-dirs=%s' % ' '.join(lapack_libs.directories),
            '--with-blas-libs=%s' % ' '.join(blas_libs.names),
            '--with-blas-lib-dirs=%s' % ' '.join(blas_libs.directories)
        ]

        if '+mpi' in self.spec:
            configure_args.append('CC=%s' % spec['mpi'].mpicc)
            configure_args.append('CXX=%s' % spec['mpi'].mpicxx)
            configure_args.append('F77=%s' % spec['mpi'].mpif77)
            configure_args.append('--with-MPI')
        else:
            configure_args.append('CC=%s' % spack_cc)
            configure_args.append('CXX=%s' % spack_cxx)
            configure_args.append('F77=%s' % spack_f77)
            configure_args.append('--without-MPI')

        cflags = ' '.join(spec.compiler_flags['cflags'])
        cxxflags = ' '.join(spec.compiler_flags['cxxflags'])
        fflags = ' '.join(spec.compiler_flags['fflags'])

        if '+openmp' in self.spec:
            configure_args.append('--with-openmp')
            cflags += ' ' + self.compiler.openmp_flag
            configure_args.append('LDFLAGS=%s' % self.compiler.openmp_flag)
        else:
            configure_args.append('--without-openmp')
            
        configure_args.append('CFLAGS=%s' % cflags)
        configure_args.append('CXXFLAGS=%s' % cxxflags)
        configure_args.append('FCFLAGS=%s' % fflags)

        if '+int64' in self.spec:
            configure_args.append('--enable-bigint')
        else:
            configure_args.append('--disable-bigint')

        if '+mixedint' in self.spec:
            configure_args.append('--enable-mixedint')
        else:
            configure_args.append('--disable-mixedint')

        if '+complex' in self.spec:
            configure_args.append('--enable-complex')
        else:
            configure_args.append('--disable-complex')

        if '+shared' in self.spec:
            configure_args.append("--enable-shared")

        if '~internal-superlu' in self.spec:
            configure_args.append("--without-superlu")
            # MLI and FEI do not build without superlu on Linux
            configure_args.append("--without-mli")
            configure_args.append("--without-fei")

        if '+superlu-dist' in self.spec:
            configure_args.append('--with-dsuperlu-include=%s' %
                                  spec['superlu-dist'].prefix.include)
            configure_args.append('--with-dsuperlu-lib=%s' %
                                  spec['superlu-dist'].libs)
            configure_args.append('--with-dsuperlu')

        if '+debug' in self.spec:
            configure_args.append("--enable-debug")
        else:
            configure_args.append("--disable-debug")

        # Hypre's source is staged under ./src so we'll have to manually
        # cd into it.
        with working_dir("src"):
            configure(*configure_args)

            make()
            if self.run_tests:
                make("check")
                make("test")
                Executable(join_path('test', 'ij'))()
                sstruct = Executable(join_path('test', 'struct'))
                sstruct()
                sstruct('-in', 'test/sstruct.in.default', '-solver', '40',
                        '-rhsone')
            make("install")

    @property
    def headers(self):
        """Export the main hypre header, HYPRE.h; all other headers can be found
        in the same directory.
        Sample usage: spec['hypre'].headers.cpp_flags
        """
        hdrs = find_headers('HYPRE', self.prefix.include, recursive=False)
        return hdrs or None

    @property
    def libs(self):
        """Export the hypre library.
        Sample usage: spec['hypre'].libs.ld_flags
        """
        is_shared = '+shared' in self.spec
        libs = find_libraries('libHYPRE', root=self.prefix, shared=is_shared,
                              recursive=True)
        return libs or None
