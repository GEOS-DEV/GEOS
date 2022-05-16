from distutils.core import setup

setup(name='geosx_mesh_tools',
      version='0.1.0',
      description='Tools for managing meshes in GEOSX',
      author='Chris Sherman',
      author_email='sherman27@llnl.gov',
      packages=['geosx_mesh_tools'],
      entry_points={'console_scripts': ['convert_abaqus = geosx_mesh_tools.abaqus_converter:main',
                                        'expand_nodeset = geosx_mesh_tools.modify_mesh:main']},
      install_requires=['meshio'])
