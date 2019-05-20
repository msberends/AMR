#######################################################################
# To push new commits to the premaster branch, run:                   #
# bash git_premaster.sh "commit message"                              #
#                                                                     #
# After successful CRAN checks, merge it to the master branch with:   #
# bash git_merge.sh                                                   #
#######################################################################

echo "•••••••••••••••••••••"
echo "• Reloading package •"
echo "•••••••••••••••••••••"
Rscript -e "devtools::load_all()"
Rscript -e "devtools::document()"
Rscript -e "devtools::install(quiet = TRUE, dependencies = FALSE)"
echo
echo "•••••••••••••••••"
echo "• Building site •"
echo "•••••••••••••••••"
Rscript -e "suppressMessages(pkgdown::build_site(lazy = TRUE, examples = FALSE))"

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
git commit -a -m "website update" --quiet
# git push --quiet
git checkout master
git merge premaster
git push --quiet
git checkout premaster
