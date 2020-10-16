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


def cmake_cache_entry(name, value, comment=""):
    """Generate a string for a cmake cache variable"""

    return 'set(%s "%s" CACHE PATH "%s")\n\n' % (name, value, comment)


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

    variant('shared', default=True, description='Build Shared Libs')
    variant('caliper', default=False, description='Build Caliper support')
    variant('mkl', default=False, description='Use the Intel MKL library.')

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
    depends_on('raja@0.12.1: ~examples~exercises')
    depends_on('raja +cuda', when='+cuda')

    depends_on('umpire@4.0.1: ~examples')
    depends_on('umpire +cuda', when='+cuda')

    depends_on('chai +raja~benchmarks~examples')
    depends_on('chai +raja+cuda', when='+cuda')

    #
    # IO
    #
    depends_on('hdf5@1.10.5: +mpi', when='~vtk')

    depends_on('conduit@0.5: ~test~fortran+mpi+hdf5~hdf5_compat')

    depends_on('silo@4.10: ~fortran~silex+mpi~zlib')

    depends_on('adiak@0.2: +mpi', when='+caliper')
    depends_on('caliper@2.4: +adiak+mpi~dyninst~callpath~papi~libpfm~gotcha~sampler~sosflow', when='+caliper')

    depends_on('pugixml@1.8:')

    #
    # Math
    #
    depends_on('intel-mkl', when='+mkl')

    depends_on('parmetis@4.0.3:')
    depends_on('superlu-dist@6.3.1: +openmp+int64')

    # trilinos
    # hypre
    # suitesparse
    # petsc

    #
    # Dev tools
    #
    depends_on('uncrustify@0.67')

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
        all of the options used to configure and build Umpire.

        For more details about 'host-config' files see:
            http://software.llnl.gov/conduit/building.html

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
        cfg = open(host_config_path, "w")
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
        cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))

        # use global spack compiler flags
        cflags = ' '.join(spec.compiler_flags['cflags'])
        if cflags:
            cfg.write(cmake_cache_entry("CMAKE_C_FLAGS", cflags))

        cxxflags = ' '.join(spec.compiler_flags['cxxflags'])
        if cxxflags:
            cfg.write(cmake_cache_entry("CMAKE_CXX_FLAGS", cxxflags))

        release_flags = "-O3 -DNDEBUG"
        cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELEASE",
                                     release_flags))
        reldebinf_flags = "-O3 -g -DNDEBUG"
        cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELWITHDEBINFO",
                                     reldebinf_flags))
        debug_flags = "-O0 -g"
        cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_DEBUG", debug_flags))

        if "+cuda" in spec:
            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# Cuda\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_option("ENABLE_CUDA", True))
            cfg.write(cmake_cache_entry("CMAKE_CUDA_STANDARD", 14))

            cudatoolkitdir = spec['cuda'].prefix
            cfg.write(cmake_cache_entry("CUDA_TOOLKIT_ROOT_DIR",
                                        cudatoolkitdir))
            cudacompiler = "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc"
            cfg.write(cmake_cache_entry("CMAKE_CUDA_COMPILER", cudacompiler))

            cmake_cuda_flags = ('-restrict --expt-extended-lambda -Werror '
                                'cross-execution-space-call,reorder,'
                                'deprecated-declarations')

            archSpecifiers = ("-mtune", "-mcpu", "-march", "-qtune", "-qarch")
            for archSpecifier in archSpecifiers:
                for compilerArg in spec.compiler_flags['cxxflags']:
                    if compilerArg.startswith(archSpecifier):
                        cmake_cuda_flags += ' -Xcompiler ' + compilerArg

            if not spec.satisfies('cuda_arch=none'):
                cuda_arch = spec.variants['cuda_arch'].value
                cmake_cuda_flags += ' -arch sm_{0}'.format(cuda_arch[0])

            cfg.write(cmake_cache_string("CMAKE_CUDA_FLAGS", cmake_cuda_flags))

            cfg.write(cmake_cache_string("CMAKE_CUDA_FLAGS_RELEASE",
                                         "-O3 -Xcompiler -O3 -DNDEBUG"))
            cfg.write(cmake_cache_string("CMAKE_CUDA_FLAGS_RELWITHDEBINFO",
                                         "-O3 -g -lineinfo -Xcompiler -O3"))
            cfg.write(cmake_cache_string("CMAKE_CUDA_FLAGS_DEBUG",
                                         "-O0 -Xcompiler -O0 -g -G"))

        else:
            cfg.write(cmake_cache_option("ENABLE_CUDA", False))

        cfg.write("#{0}\n".format("-" * 80))
        cfg.write("# RAJA\n")
        cfg.write("#{0}\n\n".format("-" * 80))

        raja_dir = spec['raja'].prefix
        cfg.write(cmake_cache_entry("RAJA_DIR", raja_dir))

        cfg.write("#{0}\n".format("-" * 80))
        cfg.write("# Umpire\n")
        cfg.write("#{0}\n\n".format("-" * 80))

        if "+umpire" in spec:
            cfg.write(cmake_cache_option("ENABLE_UMPIRE", True))
            umpire_dir = spec['umpire'].prefix
            cfg.write(cmake_cache_entry("UMPIRE_DIR", umpire_dir))
        else:
            cfg.write(cmake_cache_option("ENABLE_UMPIRE", False))

        cfg.write("#{0}\n".format("-" * 80))
        cfg.write("# CHAI\n")
        cfg.write("#{0}\n\n".format("-" * 80))

        if "+chai" in spec:
            cfg.write(cmake_cache_option("ENABLE_CHAI", True))
            chai_dir = spec['chai'].prefix
            cfg.write(cmake_cache_entry("CHAI_DIR", chai_dir))
        else:
            cfg.write(cmake_cache_option("ENABLE_CHAI", False))

        cfg.write("#{0}\n".format("-" * 80))
        cfg.write("# Caliper\n")
        cfg.write("#{0}\n\n".format("-" * 80))

        if "+caliper" in spec:
            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# Caliper\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_option("ENABLE_CALIPER", True))
            caliper_dir = spec['caliper'].prefix
            cfg.write(cmake_cache_entry("CALIPER_DIR", caliper_dir))
        else:
            cfg.write(cmake_cache_option("ENABLE_CALIPER", False))

        cfg.write("#{0}\n".format("-" * 80))
        cfg.write("# Documentation\n")
        cfg.write("#{0}\n\n".format("-" * 80))
        if "+docs" in spec:
            sphinx_dir = spec['py-sphinx'].prefix
            cfg.write(cmake_cache_string('SPHINX_EXECUTABLE',
                                         os.path.join(sphinx_dir,
                                                      'bin',
                                                      'sphinx-build')))

            doxygen_dir = spec['doxygen'].prefix
            cfg.write(cmake_cache_string('DOXYGEN_EXECUTABLE',
                                         os.path.join(doxygen_dir,
                                                      'bin',
                                                      'doxygen')))

        cfg.write("#{0}\n".format("-" * 80))
        cfg.write("# addr2line\n")
        cfg.write("#{0}\n\n".format("-" * 80))
        cfg.write(cmake_cache_option('ENABLE_ADDR2LINE', '+addr2line' in spec))

        cfg.write("#{0}\n".format("-" * 80))
        cfg.write("# Other\n")
        cfg.write("#{0}\n\n".format("-" * 80))

    def cmake_args(self):
        spec = self.spec
        host_config_path = self._get_host_config_path(spec)

        options = []
        options.extend(['-C', host_config_path])

        # Shared libs
        options.append(self.define_from_variant('BUILD_SHARED_LIBS', 'shared'))

        if '~tests~examples~benchmarks' in spec:
            options.append('-DENABLE_TESTS=OFF')
        else:
            options.append('-DENABLE_TESTS=ON')

        if '~test' in spec:
            options.append('-DDISABLE_UNIT_TESTS=ON')
        elif "+tests" in spec and ('%intel' in spec or '%xl' in spec):
            warnings.warn('The LvArray unit tests take an excessive amount of'
                          ' time to build with the Intel or IBM compilers.')

        options.append(self.define_from_variant('ENABLE_EXAMPLES', 'examples'))
        options.append(self.define_from_variant('ENABLE_BENCHMARKS',
                                                'benchmarks'))
        options.append(self.define_from_variant('ENABLE_DOCS', 'docs'))

        return options
