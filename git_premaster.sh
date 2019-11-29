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
if [ -z "$3" ]; then
  # no version number set, so get latest tags to create it
  git pull --tags --quiet
  current_tag=`git describe --tags --abbrev=0 | sed 's/v//'`
  current_commit=`git describe --tags | sed 's/.*-\(.*\)-.*/\1/'`
  # combine tag (e.g. 0.1.0) and commit number (like 40) increased by 9000 to indicate beta version
  new_version="$current_tag.$((current_commit + 9000))" # results in 0.1.0.9040
  if [ -z "$new_version" ]; then
    new_version="$current_tag.9000"
    echo
    echo "** COULD NOT CREATE NEW VERSION NUMBER! **"
    echo "Are there some unpushed changes in a new tag? Then mind NEWS.md. Assuming sequence number 9000."
    echo
  fi
  # add date to 2nd line of NEWS.md when no version number was set
  sed -i -- "2s/.*/## \<small\>Last updated: $(date '+%d-%b-%Y')\<\/small\>/" NEWS.md
else
  new_version=$3
  # rmove 2nd line of NEWS.md (the last changed date)
  sed -i -- "2s/.*//" NEWS.md
  echo "Run devtools::release() or devtools::submit_cran() after this script. Non-interactive mode will not work."
fi
# set version number to DESCRIPTION and NEWS files
sed -i -- "s/^Version: .*/Version: ${new_version}/" DESCRIPTION
sed -i -- "1s/.*/# AMR ${new_version}/" NEWS.md
rm *-- || true
echo "• First 3 lines of DESCRIPTION:"
head -3 DESCRIPTION
echo
echo "• First 2 lines of NEWS.md:"
head -2 NEWS.md
echo
echo "•••••••••••••••••••••••••••••••••"
echo "• Reloading/documenting package •"
echo "•••••••••••••••••••••••••••••••••"
Rscript -e "devtools::load_all(quiet = TRUE)"
echo "• Documenting..."
Rscript -e "suppressMessages(devtools::document())"
echo
echo "••••••••••••••••••••••••••"
echo "• Updating internal data •"
echo "••••••••••••••••••••••••••"
Rscript -e "source('data-raw/internals.R')"
echo
echo "•••••••••••••••••"
echo "• Building site •"
echo "•••••••••••••••••"
echo "• Installing..."
Rscript -e "devtools::install(quiet = TRUE, dependencies = FALSE)"
Rscript -e "suppressMessages(pkgdown::build_site(lazy = $lazy, examples = FALSE))"
echo
echo "•••••••••••••••••••••••••"
echo "• List of changed files •"
echo "•••••••••••••••••••••••••"
git status --short
echo
read -p "Uploading version ${new_version}. Continue (Y/n)? " choice
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
echo "•••••••••"
echo "• Done •"
echo "••••••••"
echo

read -p "Use R-hub to simulate all CRAN checks (y/N)? " choice
case "$choice" in
  y|Y|j|J ) ;;
  * ) exit 1;;
esac
Rscript -e "rhub::check(devtools::build(), platform = c('debian-clang-devel', 'debian-gcc-devel', 'fedora-clang-devel', 'fedora-gcc-devel', 'windows-x86_64-devel', 'debian-gcc-patched', 'solaris-x86-patched', 'debian-gcc-release', 'windows-x86_64-release', 'macos-elcapitan-release', 'windows-x86_64-oldrel'))"
echo
