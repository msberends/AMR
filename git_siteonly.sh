#######################################################################
# To push new commits to the premaster branch, run:                   #
# bash git_premaster.sh "commit message"                              #
# This creates auto version numbering in DESCRIPTION and NEWS.md.     #
#                                                                     #
# After successful CRAN checks, merge it to the master branch with:   #
# bash git_merge.sh                                                   #
#                                                                     #
# To prerelease a new version number, run:                            #
# bash git_premaster.sh "v0.x.x" FALSE "0.x.x"                        #
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
