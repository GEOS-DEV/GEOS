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

class Geosx(CachedCMakePackage, CudaPackage, ROCmPackage):
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

    # phases = ['hostconfig', 'cmake', 'build', 'install']

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


    @property
    def sys_type(self):
        sys_type = str(self.spec.architecture)
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        return sys_type

    @property
    def cache_name(self):
        var = ''
        if self.spec.statisfies( '+cuda' ):
            var = '-'.join((var, 'cuda'))
        elif self.spec.statisfies( '+rocm' ):
            var = '-'.join((var, 'rocm'))
        hostname = socket.gethostname().rstrip('1234567890')
        host_config_path = "%s-%s-%s%s.cmake" % (hostname, self.sys_type(self.spec), self.spec.compiler, var)
        dest_dir = self.stage.source_path
        host_config_path = os.path.abspath(pjoin(dest_dir, host_config_path))
        return host_config_path

    def initconfig_mpi_entries(self):
        entries = super().initconfig_mpi_entries()
        entries.insert(0, cmake_cache_option('ENABLE_MPI', True))

    def initconfig_hardware_entries(self):
        entries = super().initconfig_hardware_entries() # handles generic cuda/rocm flags, doesn't handle openmp

        prefix_entries = []
        prefix_entries.append(super().define_cmake_cache_from_variant('ENABLE_OPENMP','openmp'))
        prefix_entries.append(super().define_cmake_cache_from_variant('ENABLE_CUDA','cuda'))
        prefix_entries.append(cmake_cache_string('CMAKE_CUDA_STANDARD','17'))

        archs = self.spec.variants["cuda_arch"].value
        if archs[0] != "none":
            prefix_entries.append( cmake_cache_string("CUDA_ARCH", ';'.join( f'sm_{arch}' for arch in archs ) ) )

        cmake_cuda_flags = ('-restrict --expt-extended-lambda -Werror '
                            'cross-execution-space-call,reorder,'
                            'deprecated-declarations')

        archSpecifiers = ('-mtune', '-mcpu', '-march', '-qtune', '-qarch')
        for archSpecifier in archSpecifiers:
            for compilerArg in self.spec.compiler_flags['cxxflags']:
                if compilerArg.startswith(archSpecifier):
                    cmake_cuda_flags += ' -Xcompiler ' + compilerArg

        if archs[0] != "none":
            for arch in archs:
                cmake_cuda_flags += f" -arch sm_{arch}"

        prefix_entries.append(cmake_cache_string('CMAKE_CUDA_FLAGS', cmake_cuda_flags))

        cmake_cuda_release_flags = '-O3 -Xcompiler -O3 -DNDEBUG'
        prefix_entries.append(cmake_cache_string('CMAKE_CUDA_FLAGS_RELEASE', cmake_cuda_release_flags ))
        prefix_entries.append(cmake_cache_string('CMAKE_CUDA_FLAGS_RELWITHDEBINFO', f'-g -lineinfo {cmake_cuda_release_flags}'))
        prefix_entries.append(cmake_cache_string('CMAKE_CUDA_FLAGS_DEBUG', '-g -G -O0 -Xcompiler -O0 '))

        return [ *prefix_entries, *entries ]

    def initconfig_package_entries(self):
        hdiv = "#------------------{0}".format("-" * 30)

        entries = []
        entries.append( hdiv )
        entries.append( '# Package TPLs' )
        entries.append( hdiv )

        tpls = ( ('raja', 'RAJA', True),
                 ('umpire', 'UMPIRE', True),
                 ('chai', 'CHAI', True),
                 ('hdf5', 'HDF5', True),
                 ('conduit', 'CONDUIT', True),
                 ('silo', 'SILO', self.spec.satisfies('+silo')),
                 ('adiak', 'ADIAK', self.spec.satisfies('+caliper') and self.spec['caliper'].satisfies('+adiak')),
                 ('caliper', 'CALIPER', self.spec.satisfies('+caliper')),
                 ('pugixml', 'PUGIXML', True),
                 ('vtk', 'VTK', self.spec.satisfies('+vtk')),
                 ('fmt', 'FMT', True),
                 ('metis', 'METIS', True),
                 ('parmetis', 'PARMETIS', True),
                 ('scotch', 'SCOTCH', self.spec.satisfies('+scotch')),
                 ('superlu-dist', 'SUPERLU_DIST', self.spec.satisfies('+superlu-dist')),
                 ('suite-sparse', 'SUITESPARSE', self.spec.satisfies('+suite-sparse')),
                 ('trilinos', 'TRILINOS', self.spec.satisfies('+trilinos')),
                 ('hypre', 'HYPRE', self.spec.satisfies('+hypre')),
                 ('petsc', 'PETSC', self.spec.satisfies('+petsc'))
                 )

        for tpl, cmake_name, enable in tpls:
            entries.append(cmake_cache_option( f'ENABLE_{cmake_name}', False) )
            if enable:
                entries.append(cmake_cache_path( f'{cmake_name}_DIR', self.spec[tpl].prefix) )

        if self.spec.satisfies( '^hypre+cuda' ):
            entries.append( cmake_cache_string('ENABLE_HYPRE_DEVICE', "CUDA") )
        elif self.spec.satisfies( '^hypre+rocm' ):
            entries.append( cmake_cache_string('ENABLE_HYPRE_DEVICE', "HIP") )

        if self.spec.satisfies('lai=trilinos'):
            entries.append( cmake_cache_path('GEOSX_LA_INTERFACE', 'Trilinos') )
        if self.spec.satisfies('lai=hypre'):
            entries.append( cmake_cache_path('GEOSX_LA_INTERFACE', 'Hypre') )
        if self.spec.satisfies('lai=petsc'):
            entries.append( cmake_cache_path('GEOSX_LA_INTERFACE', 'Petsc') )

        entries.append( hdiv )
        entries.append( '# BLAS/LAPACK (System) TPLs' )
        entries.append( hdiv )
        if self.spec.satisfies('^intel-oneapi-mkl'):
            entries.append(cmake_cache_option('ENABLE_MKL', True))
            entries.append(cmake_cache_path('MKL_INCLUDE_DIRS', self.spec['intel-oneapi-mkl'].prefix.include))
            entries.append(cmake_cache_string('MKL_LIBRARIES', ';'.join( self.spec['intel-oneapi--mkl'].libs ) ))
        elif self.spec.satisfies('^essl'):
            entries.append(cmake_cache_option('ENABLE_ESSL', True))
            entries.append(cmake_cache_path('ESSL_INCLUDE_DIRS', self.spec['essl'].prefix.include))
            entries.append(cmake_cache_string('ESSL_LIBRARIES', ';'.join( self.spec['essl'].libs + self.spec['cuda'].libs ) ))

            #TODO: these should be handled internall in geosx cmakelists.txt depending on the compiler
            entries.append(cmake_cache_option('FORTRAN_MANGLE_NO_UNDERSCORE', True))
            if self.spec.satisfies('+openmp'):
                entries.append(cmake_cache_string('OpenMP_Fortran_FLAGS','-qsmp=omp'))
                entries.append(cmake_cache_string('OpenMP_Fortran_LIB_NAMES',''))
        else:
            entries.append(cmake_cache_string('BLAS_LIBRARIES', ';'.join( self.spec['blas'].libs ) ))
            entries.append(cmake_cache_string('LAPACK_LIBRARIES', ';'.join( self.spec['lapack'].libs ) ))

        entries.append( hdiv )
        entries.append( '# Python' )
        entries.append( hdiv )
        entries.append( super().define_cmake_cache_from_variant( 'ENABLE_PYGEOSX','pygeosx') )
        if self.spec.satisfies('+pygeosx') :
            entries.append(cmake_cache_path('Python3_EXECUTABLE', os.path.join( self.spec['python'].prefix.bin, 'python3')))

        entries.append( hdiv )
        entries.append( '# Documentation' )
        entries.append( hdiv )

        entries.append( super().define_cmake_cache_from_variant('ENABLE_DOCS', 'doc') )
        entries.append( super().define_cmake_cache_from_variant('ENABLE_DOXYGEN', 'doc') )
        entries.append( super().define_cmake_cache_from_variant('ENABLE_SPHINX', 'doc') )

        if self.spec.satisfies( '+doc' ):
            entries.append( cmake_cache_path('SPHINX_EXECUTABLE', os.path.join(self.spec['py-sphinx'].prefix.bin, 'sphinx-build')) )
            entries.append( cmake_cache_path('DOXYGEN_EXECUTABLE', os.path.join(self.spec['doxygen'].prefix.bin, 'doxygen')) )

        entries.append( super().define_cmake_cache_from_variant('ENABLE_UNCRUSTIFY', 'dev') )
        if self.spec.satisfies( '+dev' ):
            entries.append( hdiv )
            entries.append('# Development tools\n')
            entries.append( hdiv )
            entries.append( cmake_cache_path('UNCRUSTIFY_EXECUTABLE', os.path.join(self.spec['uncrustify'].prefix.bin, 'uncrustify')) )

        # entries.append( hdiv )
        # entries.append('# addr2line\n')
        # entries.append( hdiv )
        # entries.append( super().define_cmake_cache_from_variant('ENABLE_ADDR2LINE', 'addr2line') )

        entries.append( hdiv )
        entries.append('# Other\n')
        entries.append( hdiv )

        entries.append( cmake_cache_option('ENABLE_MATHPRESSO', False) )
        entries.append( cmake_cache_option('ENABLE_XML_UPDATES', False) )

        return entries

    def cmake_args(self):
        return [ '-C', self._get_host_config_path(self.spec), '--trace', '--trace-expand' ]

