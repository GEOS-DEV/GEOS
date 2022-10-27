# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Uncrustify(Package):
    """Source Code Beautifier for C, C++, C#, ObjectiveC, Java, and others."""

    homepage = "http://uncrustify.sourceforge.net/"
    git = "https://github.com/uncrustify/uncrustify"

    version('0.70geosx', commit="401a4098bce9dcc47e024987403f2d59d9ba7bd2")

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
