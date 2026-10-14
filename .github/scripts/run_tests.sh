#!/bin/bash
# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in the Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

# Runs an R script, in a real interactive R session if AMR_TEST_INTERACTIVE is "true", and exits with its status.
# Usage: run_tests.sh <script.R>
# R CMD check and Rscript are never interactive, so code behind interactive() is otherwise never run in CI.
# In an interactive session, R reads both its commands and the answers to prompts from stdin, and it ignores -f. Therefore stdin holds one source() line, followed by an endless stream of "1", which means 'proceed' in every utils::menu() prompt of this package.
# A timeout prevents an unanswerable prompt from hanging the job.

set +e +o pipefail
script="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"

if [ "$AMR_TEST_INTERACTIVE" = "true" ]; then
  { echo "source('$script')"; yes 1; } | timeout 3600 R --interactive --no-save --quiet
  status=${PIPESTATUS[1]}
else
  timeout 3600 Rscript "$script"
  status=$?
fi

if [ "$status" = "124" ]; then
  echo "Timed out, possibly on a prompt that could not be answered"
fi
exit $status
