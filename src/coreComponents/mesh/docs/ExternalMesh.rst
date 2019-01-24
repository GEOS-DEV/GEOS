###############################################################################
Working with an external mesh
###############################################################################

Introduction
============

GEOSX provides features to run simulations on unstructured meshes.
It uses PAMELA_ to read the external meshes and it's API to write
it into the GEOSX mesh data structure.

The supported mesh format are:

- The GMSH file format (.msh).
- The MEDIT file format (.mesh)
- The ECLIPSE file format (.egrid, .grdecl)

.. _PAMELA: https://github.com/GEOSX/PAMELA
