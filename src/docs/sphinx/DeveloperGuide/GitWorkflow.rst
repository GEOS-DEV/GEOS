**************************************
GIT Workflow for Developement in GEOSX
**************************************

The GEOSX project is hosted on github `here <https://github.com/GEOSX>`__.
For instuctions on how to clone and build GEOSX, please refer to the :ref:`GettingStartedWithGEOSX`.

hey. look at this :ref:`FunctionManager`

Branches
========
The main development branch in GEOSX is named ``develop``, and all development branches should branch off of ``develop``. 
The typical naming convention for branches are:

* ``feature/branchSubName`` for branches that introduce new features,
* ``bugfix/branchSubName`` for branches that address a bugfix,
* ``docs/branchSubName`` for branches that update documentation. 

Submodules
==========
There are several submodules in the main GEOSX repository. 
When modifying code inside a submodule, you are responsible for creating a branch (preferably with the same 
name as the top level branch) within the submodule.
For instance, if you create a branch named ``feature/greatNewFeature``, and you modify files in a submodule,
you should also make a branch named ``feature/greatNewFeature`` in the submodule. 
Also, be aware that the parent to a submodule only has the hash location for that submodule, so if you make 
a commit, you should always commit the submodule first, then the parent repository to ensure that the parent 
points to the correct submodule commmit.

Pull Requests
=============
Prior to merging your changes into develop, you must create a pull request for your branch, and potentially for
any submodules that you have modified.
The github pull request allows you to select multiple reviewers, who will review the changes and provide feedback,
request changes, and finally approval for the pull request to be merged into develop.
Also as part of the Pull Request, Github sends a request to TravisCI to run integration testing on a select 
number of platforms to ensure that the changes compile successfully on supported platforms, and do not fail 
basic testing (i.e. unit tests and small integrated tests).
Once your Pull request is approved, and all tests are passed, one of the project leads will merge your pull 
request into ``develop``.

