# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2021 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

########################################################################
# `git_premaster.sh` takes 3 parameters:                               #
#   1. Commit message (character) [mandatory]                          #
#   2. Lazy website generation (logical), with TRUE only changed       #
#      files will be processed [defaults to TRUE]                      #
#   3. Version number to be used in DESCRIPTION and NEWS.md            #
#      [defaults to current tag and last commit number + 9000]         #
#                                                                      #
# To push new commits to the premaster branch, run:                    #
# bash git_premaster.sh "commit message"                               #
# This creates auto version numbering in DESCRIPTION and NEWS.md.      #
#                                                                      #
# After successful test checks, merge it to the master branch with:    #
# bash git_merge.sh                                                    #
#                                                                      #
# To prerelease a new version number, run:                             #
# bash git_premaster.sh "v1.x.x" FALSE "1.x.x"                         #
#                                                                      #
# To only update the website, run:                                     #
# bash git_siteonly.sh                                                 #
# (which is short for 'bash git_premaster.sh "website update" FALSE')  #
########################################################################

if [ -z "$1" ]; then
  echo "FATAL - no commit message"
  exit 1
fi
if [ -z "$2" ]; then
  lazy="TRUE"
else
  lazy=$2
fi

# be sure to be on premaster branch
git checkout premaster --quiet

echo "••••••••••••••••••••••••••••••••••••••••••••"
echo "• Updating package date and version number •"
echo "••••••••••••••••••••••••••••••••••••••••••••"
sed -i -- "s/^Date: .*/Date: $(date '+%Y-%m-%d')/" DESCRIPTION
if [ -z "$3" ]; then
  # no version number set, so get latest tags to create it
  git pull --tags --quiet
  current_tag=`git describe --tags --abbrev=0 | sed 's/v//'`
  if [ -z "current_tag" ]; then
    echo "FATAL - could not determine current tag"
    exit 1
  fi
  current_tag_full=`git describe --tags | sed 's/v//'`
  current_tag_dashes=`echo $current_tag_full | grep -o "[-]" | wc -l`
  if (( "$current_tag_dashes" < 1 )); then
    # so version number is like "1.0.0", commit nr is 0
    current_commit=0
    echo "---------------"
    echo "Mind NEWS.md! Assuming sequence number 9000."
    echo "---------------"
  else
    current_commit=`git describe --tags | sed 's/.*-\(.*\)-.*/\1/'`
  fi
  if [ -z "current_commit" ]; then
    echo "FATAL - could not determine last commit index number"
    exit 1
  fi
  # combine tag (e.g. 0.1.0) and commit number (like 40) increased by 9000 to indicate beta version
  new_version="$current_tag.$((current_commit + 9000))" # results in 0.1.0.9040
  # add date to 2nd line of NEWS.md when no version number was set
  sed -i -- "2s/.*/## \<small\>Last updated: $(date '+%e %B %Y')\<\/small\>/" NEWS.md
else
  # version number set in command
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
echo "R library location:" $(Rscript -e "cat(.libPaths()[1])")
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
Rscript -e "source('data-raw/_internals.R')"
echo
echo "••••••••••••••••••••"
echo "• Building package •"
echo "••••••••••••••••••••"
echo "• Building 'data-raw/AMR_latest.tar.gz'..."
Rscript -e "x <- devtools::build(path = 'data-raw', vignettes = FALSE, manual = FALSE, binary = FALSE, quiet = TRUE)"
mv data-raw/AMR_*.tar.gz data-raw/AMR_latest.tar.gz
echo "• Installing..."
Rscript -e "devtools::install(quiet = TRUE, dependencies = FALSE)"
echo
echo "•••••••••••••••••"
echo "• Building site •"
echo "•••••••••••••••••"
Rscript -e "suppressMessages(pkgdown::build_site(lazy = $lazy, examples = FALSE))"
# add the survey page
Rscript -e "source('data-raw/create_survey_page.R')"
echo
echo "•••••••••••••••••••••••••"
echo "• List of changed files •"
echo "•••••••••••••••••••••••••"
git add .
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
# save latest changes as well
git add .
# and commit
git commit -a -m "(v${new_version}) $1" --quiet
git push --quiet
echo "Comparison:"
echo "https://github.com/msberends/AMR/compare/master...premaster?view=inline"

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
Rscript -e "rhub::check(devtools::build(), platform = rhub::platforms()[!is.na(rhub::platforms()$`cran-name`), 'name'])"
echo
