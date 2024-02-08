###############################################################################
# Copyright (c) 2014-2022, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-666778
#
# All rights reserved.
#
# This file is part of Conduit.
#
# For details, see https://lc.llnl.gov/conduit/.
#
# Please also read conduit/LICENSE
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################

import glob
import os
import shutil
import socket
from os import environ as env

import llnl.util.tty as tty

from spack import *


def cmake_cache_entry(name, value, vtype=None):
    """
    Helper that creates CMake cache entry strings used in
    'host-config' files.
    """
    if vtype is None:
        if value == "ON" or value == "OFF":
            vtype = "BOOL"
        else:
            vtype = "PATH"
    return 'set({0} "{1}" CACHE {2} "")\n\n'.format(name, value, vtype)


class Magictestlib(Package):
    """Magictestlib"""

    homepage = "http://example.com/"
    url      = "http://example.com/"
    git      = "http://example.com/"

    version('1.0.0', 'c8b277080a00041cfc4f64619e31f6d6',preferred=True)

    depends_on('zlib')

    ###################################
    # build phases used by this package
    ###################################
    phases = ['hostconfig']

    def url_for_version(self, version):
        dummy_tar_path = os.path.abspath(os.path.join(__file__, "../../magictestlib.tar.gz"))
        url = "file://" + dummy_tar_path
        return url


    def _get_host_config_path(self, spec):
        sys_type = spec.architecture
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        host_config_path = "{0}-{1}-{2}-magictestlib-{3}.cmake".format(socket.gethostname(),
                                                                  sys_type,
                                                                  spec.compiler,
                                                                  spec.dag_hash())
        dest_dir = spec.prefix
        host_config_path = os.path.abspath(join_path(dest_dir,
                                                     host_config_path))
        return host_config_path

    def hostconfig(self, spec, prefix):
        """
        This method creates a mock 'host-config' file.
        """
        if not os.path.isdir(spec.prefix):
            os.mkdir(spec.prefix)

        #######################
        # Compiler Info
        #######################
        c_compiler = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]

        #######################################################################
        # Directly fetch the names of the actual compilers to create a
        # 'host config' file that works outside of the spack install env.
        #######################################################################
        # get hostconfig name
        host_cfg_fname = self._get_host_config_path(spec)

        cfg = open(host_cfg_fname, "w")
        cfg.write("##################################\n")
        cfg.write("# spack generated host-config\n")
        cfg.write("##################################\n")
        #######################
        # Spack Compiler info
        #######################
        cfg.write("#######\n")
        cfg.write("# using %s compiler spec\n" % spec.compiler)
        cfg.write("#######\n\n")
        cfg.write("# c compiler used by spack\n")
        cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", c_compiler))
        cfg.write("# cpp compiler used by spack\n")
        cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))

        ######
        # ZLIB
        cfg.write(cmake_cache_entry("ZLIB_DIR", spec['zlib'].prefix))

        #######################
        # Finish host-config
        #######################
        cfg.write("##################################\n")
        cfg.write("# end spack generated host-config\n")
        cfg.write("##################################\n")
        cfg.close()

        host_cfg_fname = os.path.abspath(host_cfg_fname)
        tty.info("spack generated host-config file: " + host_cfg_fname)

