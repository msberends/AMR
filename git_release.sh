#######################################################################
# To push new commits to the premaster branch, run:                   #
# bash git_premaster.sh "commit message"                              #
#                                                                     #
# After successful CRAN checks, merge it to the master branch with:   #
# bash git_merge.sh                                                   #
#                                                                     #
# Initiate a new release with:                                        #
# bash git_release.sh "new_version_number"                            #
# This will edit the DESCRIPTION file and the NEWS.md file.           #
#######################################################################

if [ -z "$1" ]; then
  echo "FATAL - no version number"
  exit 1
fi

echo "••••••••••••••••••••••••••••••••••••••••••••"
echo "• Updating package date and version number •"
echo "••••••••••••••••••••••••••••••••••••••••••••"
new_version=${1}
sed -i -- "s/^Version: .*/Version: ${new_version}/" DESCRIPTION
# update 1st line of NEWS.md
sed -i -- "1s/.*/# AMR ${new_version}/" NEWS.md
echo "First 3 lines of DESCRIPTION:"
head -3 DESCRIPTION
echo
echo "First line of NEWS.md:"
head -1 NEWS.md
echo
echo "Run devtools::release() or devtools::submit_cran(). Non-interactive mode (this script) will not work."

