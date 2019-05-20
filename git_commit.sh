#######################################################################
# To push new commits to the premaster branch, run:                   #
# bash git_premaster.sh "commit message"                              #
#                                                                     #
# After successful CRAN checks, merge it to the master branch with:   #
# bash git_merge.sh                                                   #
#######################################################################

if [ -z "$1" ]; then
  echo "FATAL - no commit message"
  exit 1
fi

echo "••••••••••••••••••••••••••••••••••••••••••••"
echo "• Updating package date and version number •"
echo "••••••••••••••••••••••••••••••••••••••••••••"
sed -i -- "s/^Date: .*/Date: $(date '+%Y-%m-%d')/" DESCRIPTION
# get latest tags
git pull --tags --quiet
# get version number: latest tag + .90 + number of commits (like 0.6.1.9033)
newversion=`git describe --tags | sed 's/-/.90/' | sed 's/-.*//' | sed 's/v//'`
sed -i -- "s/^Version: .*/Version: ${newversion}/" DESCRIPTION
echo "First 3 lines of DESCRIPTION:"
head -3 DESCRIPTION
echo
echo "•••••••••••••••••••••••"
echo "• Documenting package •"
echo "•••••••••••••••••••••••"
Rscript -e "devtools::document()"
echo
echo "••••••••••••••"
echo "• Committing •"
echo "••••••••••••••"
git add .
git commit -a -m "$1" --quiet

echo
echo "Done."
