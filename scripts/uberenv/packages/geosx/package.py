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


# reset; ./scripts/uberenv/uberenv.py --spec="%gcc@8.3.1 +mkl ^chai@master"
# reset; ./scripts/uberenv/uberenv.py --spec="%clang@10.0.1 +mkl ^chai@master"

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


class Geosx(CMakePackage, CudaPackage):
    """GEOSX simulation framework."""

    homepage = "https://github.com/GEOSX/GEOSX"
    git      = "https://github.com/GEOSX/GEOSX.git"

    version('develop', branch='develop', submodules='True')

    variant('shared', default=True, description='Build Shared Libs.')
    variant('caliper', default=False, description='Build Caliper support.')
    variant('mkl', default=False, description='Use the Intel MKL library.')
    variant('suite-sparse', default=True, description='Build SuiteSparse support.')
    variant('trilinos', default=True, description='Build Trilinos support.')
    variant('hypre', default=True, description='Build HYPRE support.')
    variant('petsc', default=True, description='Build PETSc support.')
    variant('lai', default='trilinos', description='Linear algebra interface.',
            values=('trilinos', 'hypre', 'petsc'), multi=False)
    # variant('tests', default=True, description='Build tests')
    # variant('benchmarks', default=False, description='Build benchmarks')
    # variant('examples', default=False, description='Build examples')
    # variant('docs', default=False, description='Build docs')
    # variant('addr2line', default=True,
    #         description='Build support for addr2line.')

    depends_on('cmake@3.8:', type='build')
    depends_on('cmake@3.9:', when='+cuda', type='build')

    #
    # Performance portability
    #
    depends_on('raja@0.12.1: +openmp +shared ~examples ~exercises')
    depends_on('raja +cuda', when='+cuda')

    depends_on('umpire@4.0.1: ~c +shared +openmp ~examples')
    depends_on('umpire +cuda', when='+cuda')

    depends_on('chai +shared +raja ~benchmarks ~examples')
    depends_on('chai +cuda', when='+cuda')

    #
    # IO
    #
    depends_on('hdf5@1.10.5: +shared +pic +mpi', when='~vtk')

    depends_on('conduit@0.5: +shared ~test ~fortran +mpi +hdf5 ~hdf5_compat')

    depends_on('silo@4.10: ~fortran +shared ~silex +pic +mpi ~zlib')

    depends_on('adiak@0.2: +mpi +shared', when='+caliper')
    depends_on('caliper@2.4: +shared +adiak +mpi ~callpath ~papi ~libpfm ~gotcha ~sampler', when='+caliper')

    depends_on('pugixml@1.8: +shared')

    #
    # Math
    #
    depends_on('intel-mkl +shared ~ilp64', when='+mkl')

    depends_on('parmetis@4.0.3: +shared +int64')
    depends_on('superlu-dist@6.3.1: +int64 +openmp +shared')

    depends_on('suite-sparse@5.7.2: +pic +openmp', when='+suite-sparse')

    trilinos_build_options = '~fortran +openmp +shared'
    trilinos_tpls = '~boost ~glm ~gtest ~hdf5 ~hypre ~matio ~metis +mpi ~mumps ~netcdf ~suite-sparse'
    trilinos_packages = '+amesos +aztec +epetra +epetraext +ifpack +kokkos +ml +stk +stratimikos +teuchos +tpetra ~amesos2 ~anasazi ~belos ~exodus ~ifpack2 ~muelu ~sacado ~zoltan ~zoltan2'
    depends_on('trilinos@12.18.1: ' + trilinos_build_options + trilinos_tpls + trilinos_packages, when='+trilinos')

    depends_on('hypre@2.20.0: +shared +superlu-dist +mixedint +mpi +openmp', when='+hypre')
 
    petsc_build_options = '+shared +mpi'
    petsc_tpls = '+metis ~hdf5 ~hypre ~superlu-dist +int64'
    depends_on('petsc@3.13.0: ' + petsc_build_options + petsc_tpls, when='+petsc')

    #
    # Dev tools
    #
    depends_on('uncrustify@0.71:')

    #
    # Documentation
    #
    depends_on('doxygen@1.8.13:', when='+docs', type='build')
    depends_on('py-sphinx@1.6.3:', when='+docs', type='build')

    phases = ['hostconfig', 'cmake', 'build', 'install']

    @run_after('build')
    @on_package_attributes(run_tests=True)
    def check(self):
        """Searches the CMake-generated Makefile for the target ``test``
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
        host_config_path = "%s-%s-%s%s.cmake" % (hostname,
                                                 self._get_sys_type(spec),
                                                 spec.compiler, var)

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

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# OpenMP\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_option('ENABLE_OPENMP', True))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Cuda\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+cuda' in spec:
                cfg.write(cmake_cache_option('ENABLE_CUDA', True))
                cfg.write(cmake_cache_entry('CMAKE_CUDA_STANDARD', 14))

                cudatoolkitdir = spec['cuda'].prefix
                cfg.write(cmake_cache_entry('CUDA_TOOLKIT_ROOT_DIR',
                                            cudatoolkitdir))
                cudacompiler = '${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc'
                cfg.write(cmake_cache_entry('CMAKE_CUDA_COMPILER', cudacompiler))

                cmake_cuda_flags = ('-restrict --expt-extended-lambda -Werror '
                                    'cross-execution-space-call,reorder,'
                                    'deprecated-declarations')

                archSpecifiers = ('-mtune', '-mcpu', '-march', '-qtune', '-qarch')
                for archSpecifier in archSpecifiers:
                    for compilerArg in spec.compiler_flags['cxxflags']:
                        if compilerArg.startswith(archSpecifier):
                            cmake_cuda_flags += ' -Xcompiler ' + compilerArg

                if not spec.satisfies('cuda_arch=none'):
                    cuda_arch = spec.variants['cuda_arch'].value
                    cmake_cuda_flags += ' -arch sm_{0}'.format(cuda_arch[0])

                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS', cmake_cuda_flags))

                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELEASE',
                                            '-O3 -Xcompiler -O3 -DNDEBUG'))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELWITHDEBINFO',
                                            '-O3 -g -lineinfo -Xcompiler -O3'))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_DEBUG',
                                            '-O0 -Xcompiler -O0 -g -G'))

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

            io_tpls = (('hdf5', 'HDF5', True),
                    ('conduit', 'CONDUIT', True),
                    ('silo', 'SILO', True),
                    ('adiak', 'ADIAK', '+caliper' in spec),
                    ('caliper', 'CALIPER', '+caliper' in spec),
                    ('pugixml', 'PUGIXML', True))
            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# IO TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            for tpl, cmake_name, enable in io_tpls:
                if enable:
                    cfg.write(cmake_cache_entry('{}_DIR'.format(cmake_name), spec[tpl].prefix))
                else:
                    cfg.write(cmake_cache_option('ENABLE_{}'.format(cmake_name), False))

            math_tpls = (('metis', 'METIS', True),
                        ('parmetis', 'PARMETIS', True),
                        ('superlu-dist', 'SUPERLU_DIST', True),
                        ('suite-sparse', 'SUITESPARSE', '+suite-sparse' in spec),
                        ('trilinos', 'TRILINOS', '+trilinos' in spec),
                        ('hypre', 'HYPRE', '+hypre' in spec),
                        ('petsc', 'PETSC', '+petsc' in spec))
            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Math TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            for tpl, cmake_name, enable in math_tpls:
                if enable:
                    cfg.write(cmake_cache_entry('{}_DIR'.format(cmake_name), spec[tpl].prefix))
                else:
                    cfg.write(cmake_cache_option('ENABLE_{}'.format(cmake_name), False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# System Math Libraries\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+mkl' in spec:
                cfg.write(cmake_cache_option('ENABLE_MKL', True))
                cfg.write(cmake_cache_entry('MKL_INCLUDE_DIRS', spec['intel-mkl'].prefix.include))
                cfg.write(cmake_cache_list('MKL_LIBRARIES', spec['intel-mkl'].libs))


            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Documentation\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+docs' in spec:
                sphinx_dir = spec['py-sphinx'].prefix
                cfg.write(cmake_cache_entry('SPHINX_EXECUTABLE',
                                            os.path.join(sphinx_dir,
                                                        'bin',
                                                        'sphinx-build')))

                doxygen_dir = spec['doxygen'].prefix
                cfg.write(cmake_cache_entry('DOXYGEN_EXECUTABLE',
                                            os.path.join(doxygen_dir,
                                                        'bin',
                                                        'doxygen')))
            else:
                cfg.write(cmake_cache_option('ENABLE_DOCS', False))
            
            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Development tools\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+uncrustify' in spec:
                cfg.write(cmake_cache_entry('UNCRUSTIFY_EXECUTABLE',
                                            os.path.join(spec['uncrustify'].prefix, 'bin', 'uncrustify')))
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
        pass
        # spec = self.spec
        # host_config_path = self._get_host_config_path(spec)

        # options = []
        # options.extend(['-C', host_config_path])

        # # Shared libs
        # options.append(self.define_from_variant('BUILD_SHARED_LIBS', 'shared'))

        # if '~tests~examples~benchmarks' in spec:
        #     options.append('-DENABLE_TESTS=OFF')
        # else:
        #     options.append('-DENABLE_TESTS=ON')

        # if '~test' in spec:
        #     options.append('-DDISABLE_UNIT_TESTS=ON')
        # elif "+tests" in spec and ('%intel' in spec or '%xl' in spec):
        #     warnings.warn('The LvArray unit tests take an excessive amount of'
        #                   ' time to build with the Intel or IBM compilers.')

        # options.append(self.define_from_variant('ENABLE_EXAMPLES', 'examples'))
        # options.append(self.define_from_variant('ENABLE_BENCHMARKS',
        #                                         'benchmarks'))
        # options.append(self.define_from_variant('ENABLE_DOCS', 'docs'))

        # return options
