.. _GitWorkflow:

**************************************
Git Workflow
**************************************

The GEOSX project is hosted on github `here <https://github.com/GEOSX>`__.
For instructions on how to clone and build GEOSX, please refer to the :ref:`GettingStartedWithGEOSX`.
Consider consulting `https://try.github.io/ <https://try.github.io/>`_ for practical references on how to use git.

Git Credentials
=======================================

Those who want to contribute to GEOSX should setup SSH keys for authentication, and connect
to github through SSH as discussed `in this article <https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh>`_.
Before going further, you should `test your ssh connection <https://help.github.com/en/github/authenticating-to-github/testing-your-ssh-connection>`_.
If it fails (perhaps because of your institution's proxy),
you may consider the `personal access token option <https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line>`_ as an alternative.

Downloading The Code
=======================================

Once you have created an ``ssh-key`` and you have added it to your `Github` account you can download
the code through SSH.

.. code-block:: sh

   git clone git@github.com:GEOSX/GEOSX.git your_geosx_dir <--( last argument is the directory to 
                                                                place the repository in)
   cd your_geosx_dir
   git lfs install
   git submodule init
   git submodule update

If all goes well, you should have a complete copy of the GEOSX source at this point.
The most common errors people encounter here have to do with Github not recognizing
thier authentication settings.

Branching Model
===============
The current branching model used in GEOSX (starting with version 0.2) is essentially the 
`Oneflow branching model and workflow <https://www.endoflineblog.com/oneflow-a-git-branching-model-and-workflow#when-not-to-use-oneflow/>`_
with some modifications to the merging strategy, and the treatment of release 
branches. 
In GEOSX, the ``develop`` branch serves as the main branch in the ``Oneflow`` model.
The remaining branch types are described in the following subsections.

.. note::
   The early commits in GEOSX (up to version 0.2) used a 
   `Gitflow <https://nvie.com/posts/a-successful-git-branching-model/>`_
   approach for merging feature branches into develop.
   This was done without cleaning the commit history in each feature
   branch prior to the merge into develop, resulting in an overly verbose 
   (i.e. useless) history.
   Furthermore, as would be expected, having many active feature branches resulted
   in a fairly wide/spaghetti history.
   As some point in the development process, we choose try switch primarily to a 
   squash merge approach which results in a linear develop history.
   While this is fixes the spaghetti history, we do potentially lose important
   commit history during the development process.
   Options for merging are discussed in the following sections.

.. _Feature_Branches:

Feature Branches
----------------
New developments (i.e. new features, or modifications to features) are branched off 
of ``develop`` into a ``feature`` branch.
The naming of feature branches should follow ``feature/[developer]/[branch-description]``
if you expect that only a single developer will contribute to the branch, 
or ``feature/[branch-description]`` if you expect it will be a collaborative effort.
For example, if a developer named ``neo`` were to add or modify a code feature
expecting that he/she would be the only contributor, he/she would create a branch 
using the following commands to create the local branch and push it to the remote
repository:

.. code-block:: sh

  git checkout -b feature/neo/freeYourMind
  git push -u origin feature/neo/freeYourMind

However if the branch is a collaborative branch amongst many developers, the 
appropriate commands would be:

.. code-block:: sh

  git checkout -b feature/freeYourMind
  git push -u origin feature/freeYourMind
  
When ``feature`` branches are ready to be merged into ``develop``, a ``Pull Request``
should be created to perform the review and merging process. 
See below for details about :ref:`Submitting_a_Pull_Request`.

Bugfix Branches
---------------
Bugfix branches are used to fix bugs that are present in the ``develop`` branch.
A similar naming convention to that of the ``feature`` branches are used, replacing
"feature" with "bugfix" (i.e. ``bugfix/neo/squashAgentSmith``).
Typically bugfix branches are completed by a single contributor, but just as with
the ``feature`` branches, a collaborative effort may be required resulting a 
dropping the developer name from the branch name.

When ``bugfix`` branches are ready to be merged into ``develop``, a ``Pull Request``
should be created to perform the review and merging process. 
See below for details about :ref:`Submitting_a_Pull_Request`.


Release Branches
----------------
When ``develop`` has progressed to a point where we would like to create a new 
``release``, we will create a ``release candidate`` branch with the name consisting 
of ``major.minor.x`` number, where the ``x`` represents the sequence of patch tags that
will be applied to the branch.
For instance if we were releasing version ``1.2.0``, we would name the branch
``release/1.2.x``.
In contrast to the ``OneFlow`` approach, we do not delete the ``release`` branch.
Once the release candidate release is ready, it is tagged and merged back into ``develop``. 
From that point the ``release`` branch exists to provide a basis for maintaining 
the major/minor release via the patching index.
Once the ``release/major.minor.x`` branch is tagged, it only exists as a starting point 
for patches. 
Each patch on the branch is tagged using an increment in the patch number, 
resulting in a sequence of tags like ``release/1.2.0, release/1.2.1, release/1.2.2``, etc. 

.. note::
   In our git workflow, we do not close the release loop by merging the final 
   release commit into ``develop``. 
   Once the initial ``release/x.y.z`` branch is tagged ``release/x.y.0`` and merged 
   into ``develop``, no further efforts are made to merge subsequent modifications 
   in the ``release`` branch into ``develop``.

Hotfix Branches
---------------
A ``hotfix`` branch fixes a bug in the ``release`` branch.
It uses the same naming convention as a ``bugfix`` branch.
The main difference with a ``bugfix`` branch is that the primary target branch is the 
``release`` branch instead of ``develop``.
As a soft policy, merging a ``hotfix`` into a ``release`` branch should result in 
a patch increment for the release sequence of tags.
So if a ``hotfix`` was merged into ``release/1.2.z`` with a most recent tag of
``1.2.1``, the merged commit would be tagged with ``1.2.2``.
In addition to merging into the ``release`` branch, the ``hotfix`` may also target 
``develop`` if it appropriate.
If the bug exists in ``develop`` and there are no plans to fix the bug prior 
to the next ``major.minor`` release cycle, then the ``hotfix`` should be merged
into ``develop``.
As ``develop`` may have diverged significantly from the target ``release`` branch,
the ``hotfix`` branch may need further revisions for successful integration into
``develop``. 
As such, the procedure for integration of a hotfix would be to create a new ``bugfix`` 
branch that "cherry-picks" the ``hotfix`` commit/s used to patch ``release``, and 
apply appropriate revisions for a successful merge into ``develop``.
For example if a ``hotfix`` was contained into commit ``A``, then the following
commands would create the appropriate ``bugfix`` branch:

.. code-block:: sh

  git checkout develop
  git checkout -b bugfix/neo/freeYourMind
  git cherry-pick A
  git push -u origin bugfix/neo/freeYourMind

Additional commits to this ``bugfix`` branch may then be added prior to merging 
into ``develop``.
Some useful links for use of ``git cherry-pick`` are given here:

`Documentation for git cherry-pick <https://git-scm.com/docs/git-cherry-pick>`_

`Explanation of git cherry-pick <https://stackoverflow.com/questions/9339429/what-does-cherry-picking-a-commit-with-git-mean>`_

`Example usage of git cherry-pick <https://stackoverflow.com/questions/1670970/how-to-cherry-pick-multiple-commits>`_

Documentation Branches
----------------------
A ``docs`` branch is focused on writing/improving documentation for GEOSX.
The use of the ``docs`` branch name root applies to both sphinx documentation, 
or doxygen documentation.
The ``docs`` branch follow the same naming conventions as described in the :ref:`Feature_Branches`
section.
The html produced by a Documentation branches should be proof read using sphinx/doxygen 
prior to merging into ``develop``.


Keeping Your Branch Current
===========================
Over the course of a long development effort in a single ``feature`` branch, a 
developer may need either merge ``develop`` into their ``feature`` branch, or rebase
their ``feature`` branch on ``develop``.
We don't have a mandate on how you keep your branch current, but we do have 
guidelines on the branch history when merging your branch into ``develop``.
Typically, merging ``develop`` into your branch is the easiest approach, but will
lead to a complex relationship with ``develop`` with multiple interactions...which 
can lead to a confusing history.
Conversely, rebasing your branch onto ``develop`` is more difficult, but will lead 
to a linear history within the branch.
For a complex history, we will perform a squash merge into ``develop``, thereby 
the work from the branch will appear as a single commit in ``develop``.
For clean branch histories where the individual commits are meaningful and should 
be preserved, we have to option to perform a merge commit in with the PR is merged 
into ``develop``, with the addition of a merge commit, thus maintaining the commit history.


Branching Off of a Branch
===========================
During the development processes, sometimes it is appropriate to create a branch
off of a branch. 
For instance, if there is a large collaborative development effort on the branch 
``feature/theMatrix``, and a developer would like to add a self-contained and easily 
reviewable contribution to that effort, he/she should create a branch as follows:

.. code-block:: sh

  git checkout feature/theMatrix
  git checkout -b feature/smith/dodgeBullets
  git push -u origin feature/smith/dodgeBullets

If ``feature/smith/dodgeBullets`` is intended to be merged into ``feature/theMatrix``,
and the commit history of ``feature/theMatrix`` is not changed via ``git rebase``, then
the process of merging the changes back into ``feature/theMatrix`` is fairly standard.

However, if ``feature/theMatrix`` is merged into ``develop`` via a ``squash merge``, 
and then ``smith`` would like to merge ``feature/smith/dodgeBullets`` into ``develop``,
there is a substantial problem due to the diverged history of the branches.
Specifically, ``feature/smith/dodgeBullets`` branched off a commit in ``feature/theMatrix``
that doesn't exist in ``develop`` (because it was squash merged).
For simplicity, lets assume the commit hash that ``feature/smith/dodgeBullets`` 
originated from is ``CC``, and that there were commits ``CA, CB, CC, CD`` in ``feature/theMatrix``.
When ``feature/theMatrix`` was squash-merged, all of the changes appear in ``develop`` as commit ``G``.
To further complicate the situation, perhaps a complex PR was merged after ``G``, resulting
in ``E`` on develop.
The situation is illustrated by:

.. code-block:: sh

   A----B----C----D----E----F----G----E (develop)
              \                 /
               CA---CB---CC---CD        (feature/theMatrix)
                          \
                          CCA--CCB--CCC (feature/smith/dodgeBullets)

In order to successfully merge ``feature/smith/dodgeBullets`` into ``develop``, all
commits present in ``feature/smith/dodgeBullets`` after ``CC`` must be included, while discarding
``CA, CB``, which exist in ``feature/smith/dodgeBullets`` as part of its history, but not 
in ``develop``.

One "solution" is to perform a ``git rebase --onto`` of ``feature/smith/dodgeBullets`` onto
``develop``. 
Specifically, we would like to rebase ``CCA, CCB, CCC`` onto `G`, and proceed with our
development of ``feature/smith/dodgeBullets``.
This would look like:

.. code-block:: sh

   git checkout develop
   git pull
   git checkout feature/smith/dodgeBullets
   git rebase -onto G CC

As should be apparent, we have specified the starting point as ``G``, and the point 
at which we replay the commits in ``feature/smith/dodgeBullets`` as all commits 
AFTER ``CC``.
The result is:

.. code-block:: sh

   A----B----C----D----E----F----G----E (develop)
                                  \
                                 CCA'--CCB'--CCC' (feature/smith/dodgeBullets)

Now you may proceed with standard methods for keeping ``feature/smith/dodgeBullets`` 
current with ``develop``.

.. _Submitting_a_Pull_Request:

Submitting a Pull Request
======================================
Once you have created your branch and pushed to Github, you can now create a 
`Pull Request <https://github.com/GEOSX/GEOSX/pulls>`_ on Github.
The PR creates a central place to review and discuss the ongoing work on the branch. 
Creating a pull request early in the development process is preferred as it allows 
for developers to collaborate on the branch more readily.

.. note::
   When initially creating a pull request (PR) on GitHub, create it as a draft PR while
   work is ongoing and the PR isn't ready for testing, review, and merge consideration.

When you create the initial draft PR, please ensure that you apply appropriate labels.
Applying labels allows other developers to more quickly filter the live PRs and access
those that are relevant to them. Always add the `new` label upon PR creation, as well
as to the appropriate `type`, `priority`, and  `effort` labels. In addition, please
also add any appropriate `flags`.


.. note::
   If your branch and PR will resolve any open issues, be sure to `link` them to
   the PR to ensure they are appropriate resolved once the PR is merged.
   In order to `link` the issue to the PR for automatic resolution, you must use 
   one of the keywords followed by the issue number (e.g. resolves #1020) in either
   the main description of the PR, or a commit message. 
   Entries in PR comments that are not the main description or a commit message 
   will be ignored, and the issue will not be automatically closed.
   A complete list of keywords are:

   - close
   - closes
   - closed
   - fix
   - fixes
   - fixed
   - resolve
   - resolves
   - resolved
   
   For more details, see the `Github Documentation <https://docs.github.com/en/github/managing-your-work-on-github/linking-a-pull-request-to-an-issue#linking-a-pull-request-to-an-issue-using-a-keyword>`_.

Once you are satisfied with your work on the branch, you may promote the PR out of
draft status, which will allow our integrated testing suite to execute on the PR branch
to ensure all tests are passing prior to merging.

Once the tests are passing -- or in some cases immediately -- add the `flag: ready for review`
label to the PR, and be sure to tag any relevant developers to review the PR. The PR
*must* be approved by reviewers in order to be merged.

Note that whenever a pull request is merged into ``develop``, commits are either
``squashed``, or preserved depending on the cleanliness of the history.


Keeping Submodules Current
=======================================
Whenever you switch between branches locally, pull changes from ``origin`` and/or
``merge`` from the relevant branches, it is important to update the submodules to
move the ``head`` to the proper ``commit``.

.. code-block:: sh

  git submodule update --recursive

You may also wish to modify your `git pull` behavior to update your submodules
recursively for you in one command, though you forfeit some control granularity
to do so. The method for accomplishing this varies between git versions, but
as of git 2.15 you should be able to globally configure git to accomplish this via:

.. code-block:: sh

   git config --global submodule.recurse true

In some cases, code changes will require to rebaseline the ``Integrated Tests``.
If that is the case, you will need to modify the ``integrated tests submodule``.
Instructions on how to modify a submodule are presented in the following section.

Working on the Submodules
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
