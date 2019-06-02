#######################################################################
# To push new commits to the premaster branch, run:                   #
# bash git_premaster.sh "commit message"                              #
#                                                                     #
# After successful CRAN checks, merge it to the master branch with:   #
# bash git_merge.sh                                                   #
#######################################################################

# stash current changes
# git stash --quiet

# go to master
git add .
git commit -a -m "website update" --quiet
git checkout master --quiet
echo "changed branch to master"
# import everything from premaster
git merge premaster --quiet
# and send it to git
git push --quiet
echo "pushed changes to master"
# return to premaster
git checkout premaster --quiet
echo "changed branch back to premaster"
git status --short


# and get stashed changes back
# git stash apply --quiet

