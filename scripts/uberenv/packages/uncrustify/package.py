# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Uncrustify(Package):
    """Source Code Beautifier for C, C++, C#, ObjectiveC, Java, and others."""

    homepage = "https://github.com/uncrustify/uncrustify"
    git      = "https://github.com/uncrustify/uncrustify.git"

    version('master', branch='master')
    version('0.71.0', tag='uncrustify-0.71.0')
    version('0.70.1', tag='uncrustify-0.70.1')
    version('0.70.0', tag='uncrustify-0.70')
    version('0.69.0', tag='uncrustify-0.69.0')
    version('0.68.1', tag='uncrustify-0.68.1')
    version('0.68',   tag='uncrustify-0.68')
    version('0.67',   tag='uncrustify-0.67')
    version('0.66.1', tag='uncrustify-0.66.1')
    version('0.66',   tag='uncrustify-0.66')
    version('0.65',   tag='uncrustify-0.65')
    version('0.64',   tag='uncrustify-0.64')
    version('0.63',   tag='uncrustify-0.63')
    version('0.62',   tag='uncrustify-0.62')
    version('0.61',   tag='uncrustify-0.61')

    depends_on('cmake', type='build', when='@0.64:')

    @when('@0.64:')
    def install(self, spec, prefix):
        with working_dir('spack-build', create=True):
            cmake('..', *std_cmake_args)
            make()
            make('install')

    @when('@:0.62')
    def install(self, spec, prefix):
        configure('--prefix={0}'.format(self.prefix))
        make()
        make('install')
