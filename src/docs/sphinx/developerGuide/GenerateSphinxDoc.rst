Generate Sphinx documentation files
--------------------------------------------------------------------------------

- To generate the documentation files, you will need to install Sphinx using

  .. code-block:: sh

   sudo apt install python-sphinx

  Then you can generate the documentation files with the following command

 .. code-block:: sh

  cd GEOSX/build-your-platform-release
  make geosx_docs

- That will create a new folder

  .. code-block:: sh

   GEOSX/build-your-platform-release/html/docs/sphinx

which contains all the html files generated.
