.. _contribute:

=================
How to contribute
=================

Discussion of ASE development takes place on the
:ref:`ase-user <contact>` mailing list, the ``#ase``
:ref:`Matrix channel <contact>`, and on Gitlab_.

We welcome new developers who would like to help work on improving
ASE.  If you would like to contribute, you should first tell us what
you want to work on.  Use the :ref:`mailing list <contact>` for that.


GitLab repository
=================

All work on the source code takes place on Gitlab_.

.. _Gitlab: https://gitlab.com/ase/ase


Proposed git workflow
---------------------

The workflow described here has two elements:

1. The guidelines of project git-branches as a whole.
2. The workflow of the individual developer, i.e., the standard operating
   procedure when cloning, pulling, pushing etc.

Only the latter will covered at the moment. When the developer team agree
on the former it will be added here too. One comment though:

In general:

* Never work in master branch locally or on GitLab.
* Make a new branch for whatever you are doing.  When you are done, push
  it to your own repository and make a merge request from that branch in your
  repository to official master.

The above policy ensures that the master branch means the same thing in all
repositories (official and forks).

You can learn the basics of git in several places:

* `Git Reference <https://services.github.com/>`__
* `Pro Git <https://git-scm.com/book/en/v2>`__
* `Introduction to Git with Scott Chacon of GitHub
  <https://www.youtube.com/watch?v=ZDR433b0HJY>`__
* `Tech Talk: Linus Torvalds on git
  <https://www.youtube.com/watch?v=4XpnKHJAok8>`__


Aliases
-------

These aliases are quite common::

    $ git config --global alias.st status
    $ git config --global alias.ci commit
    $ git config --global alias.co checkout
    $ git config --global alias.lol "log --pretty=oneline --abbrev-commit --graph --decorate"


The first steps as a developer
------------------------------

* Register as a user on https://about.gitlab.com/
* Find and/or change your user-name in your account setting. You will need it.

  1. Go to gitlab.com/profile
  2. Click 'Account' on the left
  3. Change your user name to be whatever you want it to be

* Follow directions on https://gitlab.com/help/ssh/README for how to generate
  and add your own public ssh-key
* Go to https://gitlab.com/ase/ase and fork the project.  The 'fork' button is
  to the right of the button with the star and to the left of the 'SSH' button.
  Forking the project give you your own personal copy of the project.

You will now have a fork situated at https://gitlab.com/your-user-name/ase/

* From here on:

  - ``upstream`` refers to git@gitlab.com:ase/ase and refers to the official
    repository  of the ase project.
  - ``origin`` refers to your copy of the ase project located at
    git@gitlab.com:your-user-name/ase or https://gitlab.com/your-user-name/ase/

  For example, ``upstream/master`` refers to the master (i.e., trunk or
  development) branch of the official ase project.

* Clone your fork to ``origin`` to your local machine::

      $ git clone git@gitlab.com:your-user-name/ase.git
      $ cd ase

* Add and track ``upstream`` in your local git (done only once for each local
  git repository)::

      $ git remote add upstream git@gitlab.com:ase/ase

  You can always check the list of remote repositories that you can obtain
  code from::

      $ git remote -v

  And you can check all available branches from the remotes that you are
  tracking::

      $ git branch -a

Making changes
--------------

Changes and/or additions can be made both directly in GitLab for small
changes (see the :ref:`small changes section<making-small-changes>`) and on a
local branch in your fork.  The preferred way using command line on a local
machine is:

1) Ensure that your master branch is in sync with upstream master (since you
   should not work in master it's ok).

  * Fetch information about heads and references from "upstream" and store it in
    local git used in later merges and checkouts::

        $ git fetch upstream

  * Switch to the local branch called 'master' that is (ideally) identical to
    the upstream master branch::

        $ git checkout master
        $ git merge upstream/master --ff-only

    If the previous command fails, then it is safe to simply reset
    your master branch to the upstream master branch with the
    ``--hard`` flag.  That will delete all local changes *and*
    extraneous commits in the current branch; so make sure (e.g.,
    ``git status``, ``git log``) that you *did* remember to check out
    the master branch *and* that you have not accidentally committed
    something here that you want to save.  And use this flag
    sparingly, as it is very powerful::

        $ git reset --hard upstream/master

    If this is first time there would be no need for hard reset, unless some time
    has passed since the cloning. Still better safe than sorry.

  * It's a good idea to keep also your own origin/master identical to
    upstream/master::

        $ git push origin master

    If this command fails, then you can try again with the ``--force`` flag.
    Same as the ``reset --hard`` git command, ``git push --force`` is powerful
    and should be used sparingly.


2) Next you can do changes and additions.

  * checkout a (new) local branch with a relevant name. Let us
    change the file contribute.rst as an example::

        $ git checkout -b add-contribute-rst

    You should typically issue this command after checking out the master
    branch
    (the new branch will be based on current *HEAD*, i.e., whatever you
    have checked out at the moment).

  * If you already have this branch from some previous work, but want to do
    new work with the same branch name then you should start by resettting it
    to current upstream/master both locally and in your GitLab account::

        $ git reset --hard upstream/master
        $ git push origin add-contribute-rst

  * Make your changes. During this stage, you should keep in mind the rule
    "Commit early and often." The next three bulleted points should be done
    many times during code editing.  Each commit should be one "unit" of work.

  * Stage the files to be committed using ``git add``::

        $ git add contribute.rst

  * Check your status::

        $ git status

  * Commit the staged changes and add commit message.  If you can summarize
    your changes succinctly, then you can use the command-line syntax::

        $ git commit -m "ENH: Add developer workflow guidelines"

    But if your changes require explanation via prose, then perhaps you should
    just execute ::

        $ git commit

    And a text editor will appear.  Please observe the following guidelines
    for writing your commit message. (stolen from
    `here <https://chris.beams.io/posts/git-commit/>`_)

    The seven rules of a great git commit message

      1. Separate subject from body with a blank line
      2. Limit the subject line to 50 characters
      3. Capitalize the subject line
      4. Do not end the subject line with a period
      5. Use the imperative mood in the subject line
      6. Wrap the body at 72 characters
      7. Use the body to explain what and why vs. how

    Read the :ref:`commit message
    section<writing-the-commit-message>` guidelines for commit messages for
    some additional ase-specific information.

  * Push commits to your GitLab repository::

        $ git push --set-upstream origin add-contribute-rst

  * Go to gitlab.com/your-user-name/ase <https://gitlab.com/your-user-name/ase>
    and click on '## branches' button (where ## is the number of branches on your
    repo)

  * Find the branch 'add-contributing-rst' and click '+ Merge Request'

  * Provide informative title and more verbose description in the
    body of the Merge Request form

  * Click the green 'Submit new merge request' button

  * For last minutes corrections that you would like to include in the
    merge request too, see :ref:`the correction
    section<Last-minute-corrections>`

  * Wait for feedback from the developer community and address concerns as
    needed by adding more commits to the 'add-contribute-rst' branch on your
    personal repository and then pushing to your gitlab repository.

  * Once the developer community is satisfied with your merge request,
    anyone with push access to gitlab.com/ase/ase <https://gitlab.com/ase/ase>
    can merge your merge request and it will now be part of the master branch

  * After the merge-request is approved, delete the branch locally::

        $ git branch -D add-contribute-rst

    and on gitlab::

        $ git push origin :add-contribute-rst
        (output)
        To git@gitlab.com:your-user-name/add-contribute-rst
        - [deleted]         add-contribute-rst


.. _Last-minute-corrections:

Adding corrections to be included in a merge request
----------------------------------------------------

If at this point you would like to make last minute corrections to your
commit (it has happened many times so don't feel too bad) then instead of
closing your own merge request and resubmit a new one you can simply
go into your branch, the one that you requested to merge the first time,
and make the changes, either directly in GitLab, see the
:ref:`small changes section<making-small-changes>`, or locally *before the
merge request has been accepted!*

Since it's the branch that is merged (not just your commit) any changes you
do to that branch will be included should the merge request be accepted::

    $ vi contribute.rst
    $ git add contribute.rst
    $ git commit
    $ git push -u origin add-contribute-rst


.. _making-small-changes:

Making small changes
--------------------

Say you want to fix a typo somewhere. GitLab has an editing feature that
can come in handy. Here are the steps to do that there:

* go to https://gitlab.com/ase/ase/
* click "Files" and find the file you want to change
* click "Edit" and fix the typo
* click "Merge Requests" and add your change from the master branch
* Unless you actually want to cancel a merge request *Do NOT* click
  any buttons that reads 'Close'!

At this point someone will take a look at your change and merge it to the
official repository if the change looks good.


.. _writing-the-commit-message:

Writing the commit message
--------------------------

Commit messages should be clear and follow a few basic rules.  Example::

   ENH: add functionality X to ase.<submodule>

   The first line of the commit message starts with a capitalized
   acronym (options listed below) indicating what type of commit this
   is.  Then a blank line, then more text if needed.  Lines shouldn't be
   longer than 72 characters. Subjects shouldn't end with a period. If
   the commit is related to a ticket, indicate that with "See #3456",
   "See ticket 3456", "Closes #3456" or similar.

Describing the motivation for a change, the nature of a bug for bug fixes or
some details on what an enhancement does are also good to include in a commit
message.  Messages should be understandable without looking at the code
changes.  A commit message like ``MAINT: fixed another one`` is an example of
what not to do; the reader has to go look for context elsewhere.

Standard acronyms to start the commit message with are:

:API: an (incompatible) API change
:BLD: change related to building ase
:BUG: bug fix
:DEP: deprecate something, or remove a deprecated object
:DEV: development tool or utility
:DOC: documentation
:ENH: enhancement
:MAINT: maintenance commit (refactoring, typos, etc.)
:REV: revert an earlier commit
:STY: style fix (whitespace, PEP8)
:TST: addition or modification of tests
:REL: related to releasing ase


Code review
===========

Before you start working on a Merge Request, *please* read our
:ref:`coding conventions`.  Please also install a linter!

Hopefully someone will look at your changes and give you some
feedback.  Maybe everything is fine and things can be merged to the official
repository right away, but there could also be some more work to do like:

* make it compatible with all supported Pythons (see
  :ref:`download_and_install`).
* write more comments
* fix docstrings
* write a test
* add some documentation

This code review loop is not something we have invented to prevent you from
contributing.  Such code review is practiced by virtually all software projects
that involve more than one person.  Code review should be viewed as an
opportunity for you to learn how to write code that fits into the ASE codebase.
