#######################################################################
# To push new commits to the premaster branch, run:                   #
# bash git_premaster.sh "commit message"                              #
#                                                                     #
# After successful CRAN checks, merge it to the master branch with:   #
# bash git_merge.sh                                                   #
#######################################################################

bash git_premaster.sh "website update" FALSE

echo
echo "••••••••••••••••••••••••••••••"
echo "• Uploading to master branch •"
echo "••••••••••••••••••••••••••••••"
git checkout master
git merge premaster
git push --quiet
git checkout premaster
