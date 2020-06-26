.. _GitWorkflow:

**************************************
Git Workflow
**************************************

Git credentials
=======================================

Those who want to contribute to geosx should setup SSH keys for authentication, and connect
to github through SSH as discussed`in this article <https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh>`_.
Before going further, you should `test your ssh connection <https://help.github.com/en/github/authenticating-to-github/testing-your-ssh-connection>`_.
If it fails (perhaps because of your institution's proxy),
you may consider the `personnal access token option <https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line>`_ as an alternative.

Downloading the code
=======================================

.. code-block:: sh

   git clone git@github.com:GEOSX/GEOSX.git
   cd GEOSX
   git lfs install
   git submodule init
   git submodule update
   cd ..

If all goes well, you should have a complete copy of the GEOSX source at this point.
The most common errors people encounter here have to do with Github not recognizing
their authentication settings.


Working on the main code
=======================================

.. code-block:: sh

  git --version

In some cases we will need to rebaseline the ``Integrated Tests`` after having performed some code changes in
the ``geosx repository``. If that is the case, we will need to modify the ``integrated tests submodule``.
Instructions on how to modify a submodule are presented in the following section.

Working on the submodules
=======================================

Sometimes it may be necessary to modify one of the submodules. In order to do so,
we need to create a pull request on the submodule repository. The following steps
can be followed in order to do so.

Move to the folder of the ``submodule`` that you intend to modify.

.. code-block:: sh

  cd submodule-folder

Currently the ``submodule`` is in detached head mode we  first move to the branch ``develop`` on the
submodule repository, pull the latest changes and then create a new branch.

.. code-block:: sh

  git checkout develop
  git pull
  git checkout -b <branch-name>

We can perform some work on this branch, `add` and `commit` the changes and then push
the newly created branch to the ``submodule repository`` on which we can eventually
create a pull request.

.. code-block:: sh

  git push --set-upstream origin <branch-name>
