# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import warnings

import socket
import os

from os import environ as env
from os.path import join as pjoin

# ./scripts/uberenv/uberenv.py --spec="%clang +mkl ^chai ^caliper+papi"

# ./scripts/uberenv/uberenv.py --spec="%clang +mkl ^raja build_type=Release ^umpire build_type=Release ^chai build_type=Release ^adiak build_type=Release ^caliper+papi build_type=Release ^pugixml build_type=Release ^parmetis build_type=Release ^superlu-dist build_type=Release ^trilinos build_type=Release"

# ./scripts/uberenv/uberenv.py --spec="%clang +essl +cuda ~petsc cuda_arch=70 ^raja build_type=Release cuda_arch=70 ^umpire build_type=Release cuda_arch=70 ^chai build_type=Release cuda_arch=70 ^adiak build_type=Release ^caliper~papi build_type=Release ^pugixml build_type=Release ^parmetis build_type=Release ^superlu-dist build_type=Release ^trilinos build_type=Release"

# PETSC doesn't compile on Lassen
# ./scripts/uberenv/uberenv.py --spec="%gcc +essl ~petsc +cuda cuda_arch=70 ^cuda@10.1.243 ^raja cuda_arch=70 ^umpire cuda_arch=70 ^chai cuda_arch=70 ^caliper~papi"


def cmake_cache_entry(name, value, comment=""):
    """Generate a string for a cmake cache variable"""

    return 'set(%s "%s" CACHE PATH "%s")\n\n' % (name, value, comment)


def cmake_cache_list(name, value, comment=""):
    """Generate a list for a cmake cache variable"""

    indent = 5 + len(name)
    join_str = '\n' + ' ' * indent
    return 'set(%s %s CACHE STRING "%s")\n\n' % (name, join_str.join(value), comment)


def cmake_cache_string(name, string, comment=""):
    """Generate a string for a cmake cache variable"""

    return 'set(%s "%s" CACHE STRING "%s")\n\n' % (name, string, comment)


def cmake_cache_option(name, boolean_value, comment=""):
    """Generate a string for a cmake configuration option"""

    value = "ON" if boolean_value else "OFF"
    return 'set(%s %s CACHE BOOL "%s")\n\n' % (name, value, comment)


class Geosx(CMakePackage, CudaPackage, ROCmPackage):
    """GEOSX simulation framework."""

    homepage = "https://github.com/GEOS-DEV/GEOS"
    git = "https://github.com/GEOS-DEV/GEOS.git"

    version('develop', branch='develop', submodules='True')

    # SPHINX_BEGIN_VARIANTS

    # variant('shared', default=True, description='Build shared Libs.')
    variant('openmp', default=True, description='Build with openmp support.')

    variant('caliper', default=True, description='Build Caliper support.')
    variant('scotch', default=False, description='Build Scotch support.')

    variant('trilinos', default=False, description='Build Trilinos support.')
    variant('hypre', default=True, description='Build HYPRE support.')
    variant('petsc', default=False, description='Build PETSc support.')
    variant('suite-sparse', default=False, description='Build SuiteSparse support.')
    variant('superlu-dist', default=False, description='Build SuperLU_dist support.')

    variant('lai',
            default='hypre',
            description='Linear algebra interface.',
            values=('trilinos', 'hypre', 'petsc'),
            multi=False)

    variant('silo', default=False, description="Build with support for silo output.")
    # variant('vtk', default=False, description="Build with support for vtk output.")

    variant('pygeosx', default=False, description='Build the GEOSX python interface.')

    # SPHINX_END_VARIANTS

    # variant('tests', default=True, description='Build tests')
    # variant('benchmarks', default=False, description='Build benchmarks')
    # variant('examples', default=False, description='Build examples')
    variant( 'dev', default=False, description='Build with optional dev tools.' )
    variant( 'docs', default=False, description='Build with doc generation support.' )


    # SPHINX_BEGIN_DEPENDS

    depends_on('cmake@3.23:', type='build')

    #
    # Virtual packages
    #
    depends_on( 'mpi' )
    depends_on( 'blas' )
    depends_on( 'lapack' )

    #
    # Immediate dependecies with only absolutely required variants/versions
    #
    depends_on( 'chai +raja' )
    depends_on( 'raja' )
    depends_on( 'umpire' )
    depends_on( 'camp' )

    depends_on( 'conduit@0.5.0: +mpi +hdf5 ~hdf5_compat' )
    depends_on( 'hdf5@1.10.5: +mpi' )
    depends_on( 'pugixml@1.8:' )
    depends_on(' fmt@10.0: cxxstd=17' )
    depends_on( 'parmetis@4.0.3: +int64' )

    #
    # Variant-specific dependencies with only required variants/versions
    #
    depends_on( 'silo@4.10: +mpi', when='+silo' )
    depends_on( 'caliper@2.4: +mpi', when='+caliper' )
    depends_on( 'scotch@6.0.9: +mpi +int64', when='+scotch' )
    depends_on( 'superlu-dist', when='+superlu-dist' )

    depends_on( 'uncrustify@0.71:', when='+dev', type='build' )
    depends_on( 'doxygen@1.8.13:', when='+docs', type='build' )
    depends_on( 'py-sphinx@1.6.3:', when='+docs', type='build' )

    #
    # Additional variants/restrictions on dependencies depending on local variants / spec restrictions.
    #
    with when( '+cuda' ):
        # we specify these to forward cuda_arch appropriately
        for cuda_arch in CudaPackage.cuda_arch_values:
            cuda_arch_variant = f'cuda_arch={cuda_arch}'
            with when( cuda_arch_variant ):
                depends_on( f'chai +cuda {cuda_arch_variant}' )
                depends_on( f'raja +cuda {cuda_arch_variant}' )
                depends_on( f'umpire +cuda {cuda_arch_variant}' )
                depends_on( f'camp +cuda {cuda_arch_variant}' )
                depends_on( f'caliper +cuda {cuda_arch_variant}' )

        depends_on( f'essl +lapackforessl +cuda', when='^essl' ) # no cuda_arch here..

    with when( '+rocm' ):
        for rocm_arch in ROCmPackage.amdgpu_targets:
            rocm_arch_variant = f'amdgpu_target={rocm_arch}'
            with when( rocm_arch_variant ):
                # we specify these to forward rocm_arch appropriately
                depends_on( f'chai +rocm  {rocm_arch_variant}' )
                depends_on( f'raja +rocm {rocm_arch_variant}' )
                depends_on( f'umpire +rocm {rocm_arch_variant}' )
                depends_on( f'camp +rocm {rocm_arch_variant}' )
                depends_on( f'caliper +rocm {rocm_arch_variant}' )

    with when( '+openmp' ):
        depends_on( 'chai +openmp' )
        depends_on( 'raja +openmp' )
        depends_on( 'umpire +openmp' )
        depends_on( 'camp +openmp' )
        depends_on( 'superlu-dist +openmp', when='+superlu-dist' )

    with when( '~openmp' ):
        # many of these default to +openmp but not all
        #  so things can get mixed and give configuration errors
        #  best just to be consistent
        depends_on( 'chai ~openmp' )
        depends_on( 'raja ~openmp' )
        depends_on( 'umpire ~openmp' )
        depends_on( 'camp ~openmp' )
        # depends_on( 'superlu-dist ~openmp', when='+superlu-dist' )

    with when( '+petsc' ):
        depends_on( 'petsc@3.13.0: +mpi' )

    with when( '+suite-sparse' ):
        depends_on('suite-sparse@5.8.1:' )

        with when( '+openmp' ):
            depends_on( 'suite-sparse +openmp' )

        with when( '+cuda' ):
            depends_on( 'suite-sparse +cuda' )

    with when( '+trilinos' ):
        depends_on( 'trilinos@12.18.1:' )

    with when( '+superlu-dist' ):
        for cuda_arch in CudaPackage.cuda_arch_values:
            cuda_arch_variant = f'cuda_arch={cuda_arch}'
            with when( cuda_arch_variant ):
                depends_on( f'superlu-dist +cuda {cuda_arch_variant}' )

        for rocm_arch in ROCmPackage.amdgpu_targets:
            rocm_arch_variant = f'amdgpu_target={rocm_arch}'
            with when( rocm_arch_variant ):
                depends_on( f'superlu-dist +rocm {rocm_arch_variant}' )

    with when( '+hypre' ):
        depends_on( 'hypre +mpi' )

        for cuda_arch in CudaPackage.cuda_arch_values:
            cuda_arch_variant = f'cuda_arch={cuda_arch}'
            with when( f'{cuda_arch_variant} ^hypre+cuda' ):
                depends_on( f'hypre +cuda +unified-memory {cuda_arch_variant}' )

        for rocm_arch in ROCmPackage.amdgpu_targets:
            rocm_arch_variant = f'amdgpu_target={rocm_arch}'
            with when( f'{rocm_arch_variant} ^hypre+rocm' ):
                depends_on( f'hypre +rocm +unified-memory {rocm_arch_variant}' )

    with when( '+pygeosx' ):
        depends_on( 'python +shared +pic' )
        depends_on( 'py-numpy@1.19: +blas +lapack' )
        depends_on( 'py-scipy@1.5.2:' )
        depends_on( 'py-mpi4py@3.0.3:' )
        depends_on( 'py-pip' )

    # SPHINX_END_DEPENDS

    #
    # Conflicts
    #
    conflicts( '+cuda +rocm' )

    conflicts('~cuda ^essl', msg='Cannot use ESSL without CUDA.')
    conflicts('~trilinos lai=trilinos', msg='To use Trilinos as the Linear Algebra Interface you must build it.')
    conflicts('~hypre lai=hypre', msg='To use HYPRE as the Linear Algebra Interface you must build it.')
    conflicts('~petsc lai=petsc', msg='To use PETSc as the Linear Algebra Interface you must build it.')

    phases = ['hostconfig', 'cmake', 'build', 'install']

    @run_after('build')
    @on_package_attributes(run_tests=True)
    def check(self):
        """
        Searches the CMake-generated Makefile for the target ``test``
        and runs it if found.
        """
        with working_dir(self.build_directory):
            ctest('-V', '--force-new-ctest-process', '-j 1')

    @run_after('build')
    def build_docs(self):
        if '+docs' in self.spec:
            with working_dir(self.build_directory):
                make('docs')

    def _get_sys_type(self, spec):
        sys_type = str(spec.architecture)
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        return sys_type

    def _get_host_config_path(self, spec):
        var = ''
        if '+cuda' in spec:
            var = '-'.join([var, 'cuda'])

        hostname = socket.gethostname().rstrip('1234567890')
        host_config_path = "%s-%s-%s%s.cmake" % (hostname, self._get_sys_type(spec), spec.compiler, var)

        dest_dir = self.stage.source_path
        host_config_path = os.path.abspath(pjoin(dest_dir, host_config_path))
        return host_config_path

    def hostconfig(self, spec, prefix, py_site_pkgs_dir=None):
        """
        This method creates a 'host-config' file that specifies
        all of the options used to configure and build GEOSX.

        Note:
          The `py_site_pkgs_dir` arg exists to allow a package that
          subclasses this package provide a specific site packages
          dir when calling this function. `py_site_pkgs_dir` should
          be an absolute path or `None`.

          This is necessary because the spack `site_packages_dir`
          var will not exist in the base class. For more details
          on this issue see: https://github.com/spack/spack/issues/6261
        """

        #######################
        # Compiler Info
        #######################
        c_compiler = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]

        #######################################################################
        # By directly fetching the names of the actual compilers we appear
        # to doing something evil here, but this is necessary to create a
        # 'host config' file that works outside of the spack install env.
        #######################################################################

        sys_type = self._get_sys_type(spec)

        ##############################################
        # Find and record what CMake is used
        ##############################################

        cmake_exe = spec['cmake'].command.path
        cmake_exe = os.path.realpath(cmake_exe)

        host_config_path = self._get_host_config_path(spec)
        with open(host_config_path, "w") as cfg:
            cfg.write("#{0}\n".format("#" * 80))
            cfg.write("# Generated host-config - Edit at own risk!\n")
            cfg.write("#{0}\n".format("#" * 80))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# SYS_TYPE: {0}\n".format(sys_type))
            cfg.write("# Compiler Spec: {0}\n".format(spec.compiler))
            cfg.write("# CMake executable path: %s\n" % cmake_exe)
            cfg.write("#{0}\n\n".format("-" * 80))

            #######################
            # Compiler Settings
            #######################

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# Compilers\n")
            cfg.write("#{0}\n\n".format("-" * 80))
            cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", c_compiler))
            cflags = ' '.join(spec.compiler_flags['cflags'])
            if cflags:
                cfg.write(cmake_cache_entry("CMAKE_C_FLAGS", cflags))

            cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))
            cxxflags = ' '.join(spec.compiler_flags['cxxflags'])
            if cxxflags:
                cfg.write(cmake_cache_entry("CMAKE_CXX_FLAGS", cxxflags))

            release_flags = "-O3 -DNDEBUG"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELEASE", release_flags))
            reldebinf_flags = "-O3 -g -DNDEBUG"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELWITHDEBINFO", reldebinf_flags))
            debug_flags = "-O0 -g"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_DEBUG", debug_flags))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# MPI\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_option('ENABLE_MPI', True))
            cfg.write(cmake_cache_entry('MPI_C_COMPILER', spec['mpi'].mpicc))
            cfg.write(cmake_cache_entry('MPI_CXX_COMPILER', spec['mpi'].mpicxx))

            if sys_type in ('linux-rhel7-ppc64le', 'linux-rhel8-ppc64le'):
                cfg.write(cmake_cache_option('ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC', True))
                if socket.gethostname().rstrip('1234567890') == "lassen":
                    cfg.write(cmake_cache_entry('MPIEXEC', 'lrun'))
                    cfg.write(cmake_cache_entry('MPIEXEC_NUMPROC_FLAG', '-n'))
                else:
                    cfg.write(cmake_cache_entry('MPIEXEC', 'jsrun'))
                    cfg.write(cmake_cache_list('MPIEXEC_NUMPROC_FLAG', ['-g1', '--bind', 'rs', '-n']))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# OpenMP\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_option('ENABLE_OPENMP', spec.satisfies('+openmp')))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Cuda\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+cuda' in spec:
                cfg.write(cmake_cache_option('ENABLE_CUDA', True))
                cfg.write(cmake_cache_entry('CMAKE_CUDA_STANDARD', 17))

                cfg.write(cmake_cache_entry('CUDA_TOOLKIT_ROOT_DIR', spec['cuda'].prefix))
                cfg.write(cmake_cache_entry('CMAKE_CUDA_HOST_COMPILER', '${CMAKE_CXX_COMPILER}'))
                cfg.write(cmake_cache_entry('CMAKE_CUDA_COMPILER', '${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc'))
                # since cuda_arch is a multi-value variant, .value is a tuple
                cfg.write(cmake_cache_string('CUDA_ARCH', ','.join( f'sm_{arch}' for arch in spec.variants['cuda_arch'].value ) ) )
                # since cuda_arch is a multi-value variant, .value is a tuple
                cfg.write(cmake_cache_string('CMAKE_CUDA_ARCHITECTURES', ','.join( spec.variants['cuda_arch'].value ) ) )

                cmake_cuda_flags = ('-restrict --expt-extended-lambda -Werror '
                                    'cross-execution-space-call,reorder,'
                                    'deprecated-declarations')

                archSpecifiers = ('-mtune', '-mcpu', '-march', '-qtune', '-qarch')
                for archSpecifier in archSpecifiers:
                    for compilerArg in spec.compiler_flags['cxxflags']:
                        if compilerArg.startswith(archSpecifier):
                            cmake_cuda_flags += ' -Xcompiler ' + compilerArg

                # since cuda_arch is a multi-value variant, .value is a tuple
                for arch in spec.variants['cuda_arch'].value:
                    cmake_cuda_flags += f" -arch sm_{arch}"

                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS', cmake_cuda_flags))

                cmake_cuda_release_flags = '-O3 -Xcompiler -O3 -DNDEBUG'
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELEASE', cmake_cuda_release_flags ))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELWITHDEBINFO', f'-g -lineinfo {cmake_cuda_release_flags}'))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_DEBUG', '-g -G -O0 -Xcompiler -O0 '))

            else:
                cfg.write(cmake_cache_option('ENABLE_CUDA', False))

            performance_portability_tpls = (('raja', 'RAJA', True),
                                            ('umpire', 'UMPIRE', True),
                                            ('chai', 'CHAI', True))
            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Performance Portability TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            for tpl, cmake_name, enable in performance_portability_tpls:
                if enable:
                    cfg.write(cmake_cache_entry('{}_DIR'.format(cmake_name), spec[tpl].prefix))
                else:
                    cfg.write(cmake_cache_option('ENABLE_{}'.format(cmake_name), False))

            have_adiak = spec.satisfies('+caliper') and '+adiak' in spec['caliper'].variants

            io_tpls = (('hdf5', 'HDF5', True),
                       ('conduit', 'CONDUIT', True),
                       ('silo', 'SILO', spec.satisfies('+silo')),
                       ('adiak', 'ADIAK', have_adiak),
                       ('caliper', 'CALIPER', spec.satisfies('+caliper')),
                       ('pugixml', 'PUGIXML', True),
                       ('vtk', 'VTK', spec.satisfies('+vtk')),
                       ('fmt', 'FMT', True))
            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# IO TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            for tpl, cmake_name, enable in io_tpls:
                if enable:
                    cfg.write(cmake_cache_entry('{}_DIR'.format(cmake_name), spec[tpl].prefix))
                else:
                    cfg.write(cmake_cache_option('ENABLE_{}'.format(cmake_name), False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# System Math Libraries\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if spec.satisfies('^intel-oneapi-mkl'):
                cfg.write(cmake_cache_option('ENABLE_MKL', True))
                cfg.write(cmake_cache_entry('MKL_INCLUDE_DIRS', spec['intel-oneapi-mkl'].prefix.include))
                cfg.write(cmake_cache_list('MKL_LIBRARIES', spec['intel-oneapi--mkl'].libs))
            elif spec.satisfies('^essl'):
                cfg.write(cmake_cache_option('ENABLE_ESSL', True))
                cfg.write(cmake_cache_entry('ESSL_INCLUDE_DIRS', spec['essl'].prefix.include))
                cfg.write(cmake_cache_list('ESSL_LIBRARIES', spec['essl'].libs + spec['cuda'].libs))

                cfg.write(cmake_cache_option('FORTRAN_MANGLE_NO_UNDERSCORE', True))
            else:
                cfg.write(cmake_cache_list('BLAS_LIBRARIES', spec['blas'].libs))
                cfg.write(cmake_cache_list('LAPACK_LIBRARIES', spec['lapack'].libs))

            # yapf: disable
            math_tpls = (
              ('metis', 'METIS', True),
              ('parmetis', 'PARMETIS', True),
              ('scotch', 'SCOTCH', spec.satisfies('+scotch')),
              ('superlu-dist', 'SUPERLU_DIST', spec.satisfies('+superlu-dist')),
              ('suite-sparse', 'SUITESPARSE', spec.satisfies('+suite-sparse')),
              ('trilinos', 'TRILINOS', spec.satisfies('+trilinos')),
              ('hypre', 'HYPRE', spec.satisfies('+hypre')),
              ('petsc', 'PETSC', spec.satisfies('+petsc'))
            )
            # yapf: enable

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Math TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            for tpl, cmake_name, enable in math_tpls:
                if enable:
                    cfg.write(cmake_cache_entry('{}_DIR'.format(cmake_name), spec[tpl].prefix))
                    if tpl == 'hypre' and spec.satisfies('^hypre+cuda'):
                        cfg.write(cmake_cache_string('ENABLE_HYPRE_DEVICE', "CUDA"))
                    elif tpl == 'hypre' and spec.satisfies('^hypre+rocm'):
                        cfg.write(cmake_cache_string('ENABLE_HYPRE_DEVICE', "HIP"))
                else:
                    cfg.write(cmake_cache_option('ENABLE_{}'.format(cmake_name), False))

            if spec.satisfies('lai=trilinos'):
                cfg.write(cmake_cache_entry('GEOSX_LA_INTERFACE', 'Trilinos'))
            if spec.satisfies('lai=hypre'):
                cfg.write(cmake_cache_entry('GEOSX_LA_INTERFACE', 'Hypre'))
            if spec.satisfies('lai=petsc'):
                cfg.write(cmake_cache_entry('GEOSX_LA_INTERFACE', 'Petsc'))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Python\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+pygeosx' in spec:
                cfg.write(cmake_cache_option('ENABLE_PYGEOSX', True))
                cfg.write(cmake_cache_entry('Python3_EXECUTABLE', os.path.join(spec['python'].prefix.bin, 'python3')))
            else:
                cfg.write(cmake_cache_option('ENABLE_PYGEOSX', False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Documentation\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if spec.satisfies('+doc'):
                sphinx_bin_dir = spec['py-sphinx'].prefix.bin
                cfg.write(cmake_cache_entry('SPHINX_EXECUTABLE', os.path.join(sphinx_bin_dir, 'sphinx-build')))

                doxygen_bin_dir = spec['doxygen'].prefix.bin
                cfg.write(cmake_cache_entry('DOXYGEN_EXECUTABLE', os.path.join(doxygen_bin_dir, 'doxygen')))
            else:
                cfg.write(cmake_cache_option('ENABLE_DOCS', False))
                cfg.write(cmake_cache_option('ENABLE_DOXYGEN', False))
                cfg.write(cmake_cache_option('ENABLE_SPHINX', False))

            if spec.satisfies('+dev'):
                cfg.write('#{0}\n'.format('-' * 80))
                cfg.write('# Development tools\n')
                cfg.write('#{0}\n\n'.format('-' * 80))
                cfg.write(cmake_cache_entry('UNCRUSTIFY_EXECUTABLE', os.path.join(spec['uncrustify'].prefix.bin, 'uncrustify')))
            else:
                cfg.write(cmake_cache_option('ENABLE_UNCRUSTIFY', False))

            # cfg.write('#{0}\n'.format('-' * 80))
            # cfg.write('# addr2line\n')
            # cfg.write('#{0}\n\n'.format('-' * 80))
            # cfg.write(cmake_cache_option('ENABLE_ADDR2LINE', '+addr2line' in spec))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Other\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            cfg.write(cmake_cache_option('ENABLE_MATHPRESSO', False))
            cfg.write(cmake_cache_option('ENABLE_XML_UPDATES', False))

    def cmake_args(self):
        return [ '-C', self._get_host_config_path(self.spec), '-