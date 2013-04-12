#!/bin/csh

# update the trunk
pushd ../../trunk/
svn update
popd

# update my branch (in case other people have made changes to it)
svn update

# merge all the recent changes made in the trunk into my branch
# (this branch was created from revision r21712, and "HEAD" indicates the current revision of the trunk
## DK DK
## DK DK  answer 'tc' in case of a detected conflict:   (tc) theirs-conflict  - accept their version for all conflicts
## DK DK
svn merge ../../trunk -r21712:HEAD

# commit the merge
svn commit -m "merged the trunk into Vadim's branch"

