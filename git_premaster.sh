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
echo "•••••••••••••••••••••••••••••••••"
echo "• Reloading/Documenting package •"
echo "•••••••••••••••••••••••••••••••••"
Rscript -e "devtools::load_all(quiet = TRUE)"
Rscript -e "suppressMessages(devtools::document())"
Rscript -e "devtools::install(quiet = TRUE, dependencies = FALSE)"
echo
echo "•••••••••••••••••"
echo "• Building site •"
echo "•••••••••••••••••"
Rscript -e "suppressMessages(pkgdown::build_site(lazy = TRUE, examples = FALSE))"

echo
echo "•••••••••••••••••••••••••"
echo "• List of changed files •"
echo "•••••••••••••••••••••••••"
git status --short
echo

read -p "Continue (Y/n)? " choice
case "$choice" in
  n|N ) exit 1;;
  * ) ;;
esac

echo
echo "•••••••••••••••••••••••••••"
echo "• Uploading to repository •"
echo "•••••••••••••••••••••••••••"
git add .
git commit -a -m "$1" --quiet
git push --quiet
echo "Comparison:"
echo "https://gitlab.com/msberends/AMR/compare/master...premaster?view=inline"

echo
echo "Done."
