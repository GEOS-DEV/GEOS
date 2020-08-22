.. _GitWorkflow:

**************************************
Git Workflow
**************************************

The GEOSX project is hosted on github `here <https://github.com/GEOSX>`__.
For instuctions on how to clone and build GEOSX, please refer to the :ref:`GettingStartedWithGEOSX`.
Consider consulting `https://try.github.io/ <https://try.github.io/>`_ for practical references on how to use git.

Git credentials
=======================================

Those who want to contribute to GEOSX should setup SSH keys for authentication, and connect
to github through SSH as discussed `in this article <https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh>`_.
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

Once you have cloned the GEOSX repository, you can start developing by creating your
own branch. GEOSX branching and merging for development is a simplification of the the
`GitFlow <https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow>`_ branching model.

In particular your branch name should follow the ``[purpose]/[developer]/[branch-description]``
naming schema, where `purpose` is usually one of ``feature, bugfix, doc``, `developer` is
your chosen identifier (often simply your github username), and branch `branch-description`
is a short description of the purpose of the branch.

Create your own branch using the following command:

.. code-block:: sh

  git checkout -b <new-branch-name>

For example, if you were going to work on a new physics solver, you might name your
branch `feature/yourName/xxx-physics-solver`.

You can now perform work on the newly created branch, ``add`` it and ``commit`` your
work. Once you have done so you can push your changes to the geosx repository and
set up an upstream branch to be tracked.

.. code-block:: sh

  git push --set-upstream origin <new-branch-name>

Once you have created your branch and pushed to remote, you can now create a pull
request on github, to create a central place to review and discuss the ongoing work
on the branch. Creating a pull request early in the development process is prefered
as it allows for developers to collaborate on the branch more readily.

Submitting a Pull Request
======================================

When initially creating a pull request (PR) on GitHub, create it is a draft PR while
work is ongoing and the PR isn't ready for testing, review, and merge consideration.

When you create the initial draft PR, please ensure that you apply appropriate labels.
Applying labels allows other developers to more quickly filter the live PRs and access
those that are relevant to them. Always add the `new` label upon PR creation, as well
as to the appropriate `type`, `priority`, and  `effort` labels. In addition, please
also add any appropriate `flags`.

Note, if your branch and PR will resolve any open issues, be sure to `link` them to
the PR to ensure they are appropriate resolved once the PR is merged.

Once you are satisfied with your work on the branch, you may promote the PR out of
draft status, which will allow our integrated testing suite to execute on the PR branch
to ensure all tests are passing prior to merging.

Once the tests are passing -- or in some cases immediately -- add the `flag: ready for review`
label to the PR, and be sure to tag any relevant developers to review the PR. The PR
*must* be approved by reviewers in order to be merged.

Note that commits are ``squashed``  whenever a pull request is merged into ``develop``.


Whenever you switch between branches locally, pull changes from ``origin`` and/or
``merge`` from the relevant branches, it is important to update the submodules to
move the ``head`` to the proper ``commit``.

.. code-block:: sh

  git submodule update --recursive

You may also wish to modify your `git pull` behavior to update your submodules
recursively for you in one command, though you forfeit some control granularity
to do so. The method for accomplishing this varies between git versions, but
as of git 2.15 you should be able to globaly configure git to accomplish this via:

.. code-block:: sh

   git config --global submodule.recurse true

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

Currently the ``submodule`` is in detached head mode, so you first need to move
to the main branch (either ``develop`` or ``master``) on the
submodule repository, pull the latest changes, and then create a new branch.

.. code-block:: sh

  git checkout <main-branch>
  git pull
  git checkout -b <branch-name>

You can perform some work on this branch, `add` and `commit` the changes and then push
the newly created branch to the ``submodule repository`` on which you can eventually
create a pull request using the same process discussed above in :ref:`Submitting a Pull Request`.

.. code-block:: sh

  git push --set-upstream origin <branch-name>


Resolving Submodule Changes In Primary Branch PRs
=================================================

When you conduct work on a submodule as described above during work on a primary GEOSX
branch which has a PR, the merging procedure requires that the submodule referenced
by the GEOSX PR branch be consistent with the submodule in the main branch of the project.
This is checked and enforced via TravisCI.

Thus in order to merge a PR that includes modifications to submodules, the various PRs for
each repository should be staged and finalized, to the point they are all ready to be merged,
with higher-level PRs in the merge hierarchy having the correct submodule references for the
current main branch for their repository.

Starting from the bottom of the submodule hierarchy, the PRs are resolved, after which the
higher-level PRs with reference to a resolved PR must update their submodule references
to point to the new main branch of the submodule with the just-resolved PR merged.
After any required automated tests pass, the higher-level PRs can then be merged.

The name of the main branch of each submodule is presented in the table below.

================    ================
Submodule           Main branch
================    ================
blt                 develop
LvArray             develop
integratedTests     develop
GEOSX_PTP           master
hdf5_interface      master
PAMELA              master
PVTPackage          master
================    ================
