# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

########################################################################
# `git_premaster.sh` takes 3 parameters:                               #
#   1. Commit message (character) [mandatory]                          #
#   2. Lazy website generation (logical), with FALSE only changed      #
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

bash git_premaster.sh "website update" FALSE

echo
echo "••••••••••••••••••••••••••••••"
echo "• Uploading to master branch •"
echo "••••••••••••••••••••••••••••••"
git checkout master
git merge premaster
git push --quiet
git checkout premaster
