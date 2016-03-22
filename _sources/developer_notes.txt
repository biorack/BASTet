Developer Notes
===============

Building the online documentation
---------------------------------

The sources for the sphinx docs are part of the doc/ folder and can be build locally as usual via, e.g., ``make html.`` This will build the documentation locally only as part of the doc/_build folder.

The online documentation is hosted via GitHub pages and are part of the gh-pages branch. The sphinx docs build scripts are conveniently set up to simply the update of the online documentation as part of the regular development process. However, before we can do this we need to do the following simple setup first:


* ``cd`` to the your local copy of the BASTet repo (i.e., where the /omsi and /doc folder are located)
* ``cd ..``
* ``mkdir bastet_docs``
* ``cd bastet_docs``
* ``git clone https://github.com/oruebel/BASTet.git html``
* ``git checkout -b gh-pages remotes/origin/ghpages``

Once we have created and setup our ``bastet_docs`` repo for the ``gh-pages`` branch we can now build the online documentation as usual from the ``doc`` folder of our development repo via the following commands:

* ``make htmlpublic`` : Rebuild the html docs locally and then copy the docs to the ../../bastet_docs repo with the gh-pages branch. This will only make local changes without commiting or publishing anything.
* ``make latexpdfpublic`` : Rebuild the latex and pdf and then copy the pdf to the ../../bastet_docs repo with the gh-pages branch.  This will only make local changes without commiting or publishing anything.
* ``make updatepublic`` : Once you have confirmed that the docs are correct, you can rebuild the html and pdf docs for publications using this command. The command also commits and pushes all changes to the docs back to GitHub for publication, so that they become immediately available online.
