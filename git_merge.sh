# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis for R                        #
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
# how to conduct AMR analysis: https://msberends.github.io/AMR/        #
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
Rscript -e "rhub::check(devtools::build(), platform = rhub::platforms()[!is.na(rhub::platforms()$`cran-name`), 'name'])"
echo

# and get stashed changes back
# git stash apply --quiet

