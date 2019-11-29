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

# stash current changes
# git stash --quiet

# go to master
git checkout master --quiet
echo "• changed branch to master"

# import everything from premaster
git merge premaster --quiet
# and send it to git
git push --quiet
echo "• pushed changes to master"

# return to premaster
git checkout premaster --quiet
echo "• changed branch back to premaster"
git status --short
echo

read -p "Use R-hub to simulate all CRAN checks (y/N)? " choice
case "$choice" in
  y|Y|j|J ) ;;
  * ) exit 1;;
esac
Rscript -e "rhub::check(devtools::build(), platform = c('debian-clang-devel', 'debian-gcc-devel', 'fedora-clang-devel', 'fedora-gcc-devel', 'windows-x86_64-devel', 'debian-gcc-patched', 'solaris-x86-patched', 'debian-gcc-release', 'windows-x86_64-release', 'macos-elcapitan-release', 'windows-x86_64-oldrel'))"
echo

# and get stashed changes back
# git stash apply --quiet

