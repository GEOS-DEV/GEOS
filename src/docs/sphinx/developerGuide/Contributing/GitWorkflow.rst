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
you may consider the `personal access token option <https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line>`_ as an alternative.

Downloading the code
=======================================

Once you have created an ``ssh-key`` and you have added it to your `Github` account you can download
the code through SSH.

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

Here, we provide only a few basic git commands that are useful to start working with
geosx. This is not at all an exhaustive guide to git and we encourage to explore the
git manual for specific doubts.

Once you have cloned the geosx repository, you can start developing by creating your
own branch.

.. code-block:: sh

  git checkout -b <new-branch-name>

You can now perform work on the newly created branch, ``add`` it and ``commit`` your
work. One you have done so you can push your changes to the geosx repository and
set up an upstream branch to be tracked.

.. code-block:: sh

  git push --set-upstream origin <new-branch-name>

You can now generate a pull request. Once you are satisfied with your work
the pull request can be reviewed and merged into the main branch. Note, that commits
are ``squashed``  whenever a pull request is merged into ``develop``.

Whenever you switch between branches, pull changes from ``origin`` or ``merge``
branches, it is important to update the submodules to move the ``head`` to the proper ``commit``.

.. code-block:: sh

  git submodule update

In some cases, code changes will require to rebaseline the ``Integrated Tests``.
If that is the case, you will need to modify the ``integrated tests submodule``.
Instructions on how to modify a submodule are presented in the following section.

Working on the submodules
=======================================

Sometimes it may be necessary to modify one of the submodules. In order to do so,
you need to create a pull request on the submodule repository. The following steps
can be followed in order to do so.

Move to the folder of the ``submodule`` that you intend to modify.

.. code-block:: sh

  cd submodule-folder

Currently the ``submodule`` is in detached head mode, so you first need to move to the branch ``develop`` on the
submodule repository, pull the latest changes, and then create a new branch.

.. code-block:: sh

  git checkout develop
  git pull
  git checkout -b <branch-name>

You can perform some work on this branch, `add` and `commit` the changes and then push
the newly created branch to the ``submodule repository`` on which you can eventually
create a pull request.

.. code-block:: sh

  git push --set-upstream origin <branch-name>
