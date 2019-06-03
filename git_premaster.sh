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
if [ -z "$2" ]; then
  lazy="TRUE"
else
  lazy=$2
fi

echo "••••••••••••••••••••••••••••••••••••••••••••"
echo "• Updating package date and version number •"
echo "••••••••••••••••••••••••••••••••••••••••••••"
sed -i -- "s/^Date: .*/Date: $(date '+%Y-%m-%d')/" DESCRIPTION
# get latest tags
git pull --tags --quiet
current_tag=`git describe --tags --abbrev=0 | sed 's/v//'`
current_commit=`git describe --tags | sed 's/.*-\(.*\)-.*/\1/'`
# combine tag (e.g. 0.1.0) and commit number (like 40) increased by 9000 to indicate beta version
new_version="$current_tag.$((current_commit + 9000))" # results in 0.1.0.9040
if [ -z "$new_version" ]; then
  new_version="$current_tag.9000"
  echo
  echo "** COULD NOT CREATE NEW VERSION NUMBER! **"
  echo "Are there some unpushed changes in a new tag?? Then mind NEWS.md. Assuming sequence number 9000."
  echo
fi
sed -i -- "s/^Version: .*/Version: ${new_version}/" DESCRIPTION
# update 1st line of NEWS.md
sed -i -- "1s/.*/# AMR ${new_version}/" NEWS.md
echo "First 3 lines of DESCRIPTION:"
head -3 DESCRIPTION
echo
echo "First line of NEWS.md:"
head -1 NEWS.md
echo
read -p "Continue (Y/n)? " choice
case "$choice" in
  n|N ) exit 1;;
  * ) ;;
esac
echo
echo "•••••••••••••••••••••••••••••••••"
echo "• Reloading/documenting package •"
echo "•••••••••••••••••••••••••••••••••"
Rscript -e "devtools::load_all(quiet = TRUE)"
Rscript -e "suppressMessages(devtools::document())"
Rscript -e "devtools::install(quiet = TRUE, dependencies = FALSE)"
echo
echo "••••••••••••••••••••••••••"
echo "• Updating internal data •"
echo "••••••••••••••••••••••••••"
Rscript -e "source('data-raw/internals.R')"
echo
echo "•••••••••••••••••"
echo "• Building site •"
echo "•••••••••••••••••"
Rscript -e "suppressMessages(pkgdown::build_site(lazy = $lazy, examples = FALSE))"
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
git commit -a -m "(v$new_version) $1" --quiet
git push --quiet
echo "Comparison:"
echo "https://gitlab.com/msberends/AMR/compare/master...premaster?view=inline"

echo
echo "Done."
