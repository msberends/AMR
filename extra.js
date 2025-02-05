/*
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
# Center Groningen in The Netherlands, in collaboration with many      #
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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #
*/

$(document).ready(function() {
  // add GPT assistant info
  $('aside').prepend('<div class="amr-gpt-assistant"><a target="_blank" href="https://chatgpt.com/g/g-M4UNLwFi5-amr-for-r-assistant"><img src="https://github.com/msberends/AMR/raw/main/pkgdown/assets/AMRforRGPT.svg"></a></div>');

  // replace 'Developers' with 'Maintainers' on the main page, and "Contributors" on the Authors page
  $(".developers h2").text("Maintainers");
  $(".citation h2:nth(0)").text("All contributors");
  $(".citation h2:nth(1)").text("How to Cite");

  // remove \donttest and \dontrun texts in Examples
  if ($("#ref-examples ~ div pre").length > 0) {
    $("#ref-examples ~ div pre").html($("#ref-examples ~ div pre").html().replaceAll("# \\donttest{", ""));
    $("#ref-examples ~ div pre").html($("#ref-examples ~ div pre").html().replaceAll("# \\dontrun{", ""));
    $("#ref-examples ~ div pre").html($("#ref-examples ~ div pre").html().replaceAll("# }", ""));
  }
  
  // rename 'Package index' to 'Manual' (weird design choice to pick 'Package index')
  if ($(".template-reference-index").length > 0) {
    $(".template-reference-index .page-header h1").text("Manual");
    document.title = document.title.replace("Package index", "Manual");
  }

  // remove leading newline in code examples on changelog
  if ($("body .template-news").length > 0) {
    $("body .template-news").html($("body .template-news").html().replaceAll('sourceCode R">\n<span', 'sourceCode R"><span'));
    $("body .template-news").html($("body .template-news").html().replaceAll('sourceCode R"><span></span>\n<span', 'sourceCode R"><span'));
  }

  // make Python part more fancy - prepare for CSS
  if (window.location.href.includes('AMR_for_Python')) {
    $('body').addClass('amr-for-python'); /* to set colours in CSS */
    $('img[src="../logo.svg"]').attr('src', '../logo_python.svg'); /* replace base logo */
    $('img[src="https://github.com/msberends/AMR/raw/main/pkgdown/assets/AMRforRGPT.svg"]').attr('src', 'https://github.com/msberends/AMR/raw/main/pkgdown/assets/AMRforRGPT_python.svg'); /* replace GPT logo */
  }
  
  // add country flags
  function country_flag(x) {
    if (typeof(x) != "undefined") {
      const prefix = '<img src="lang_';
      const suffix = '.svg" style="height: 13px !important; border: 1px solid #cccccc; aspect-ratio: 1 / 0.6; object-fit: cover; vertical-align: initial !important;"> ';
      x = x.replace("Alex", prefix + "de" + suffix + "Alex");
      x = x.replace("Andrew", prefix + "us" + suffix + "Andrew");
      x = x.replace("Anita", prefix + "au" + suffix + "Anita");
      x = x.replace("Annick", prefix + "nl" + suffix + "Annick");
      x = x.replace("Anthony", prefix + "en" + suffix + "Anthony");
      x = x.replace("Anton", prefix + "uk" + suffix + "Anton");
      x = x.replace("Bart", prefix + "nl" + suffix + "Bart");
      x = x.replace("Bhanu", prefix + "nl" + suffix + "Bhanu");
      x = x.replace("Casper", prefix + "nl" + suffix + "Casper");
      x = x.replace("Christian", prefix + "de" + suffix + "Christian");
      x = x.replace("Corinna", prefix + "nl" + suffix + "Corinna");
      x = x.replace("Dennis", prefix + "nl" + suffix + "Dennis");
      x = x.replace("Dmytro", prefix + "uk" + suffix + "Dmytro");
      x = x.replace("Eric", prefix + "nl" + suffix + "Eric");
      x = x.replace("Erwin", prefix + "nl" + suffix + "Erwin");
      x = x.replace("Gwen", prefix + "en" + suffix + "Gwen");
      x = x.replace("Jason", prefix + "ca" + suffix + "Jason");
      x = x.replace("Javier", prefix + "ca" + suffix + "Javier");
      x = x.replace("Jonas", prefix + "de" + suffix + "Jonas");
      x = x.replace("Judith", prefix + "nl" + suffix + "Judith");
      x = x.replace("Larisse", prefix + "sa" + suffix + "Larisse");
      x = x.replace("Matthew", prefix + "ca" + suffix + "Matthew");
      x = x.replace("Matthijs", prefix + "nl" + suffix + "Matthijs");
      x = x.replace("Peter", prefix + "en" + suffix + "Peter");
      x = x.replace("Rogier", prefix + "nl" + suffix + "Rogier");
      x = x.replace("Sofia", prefix + "sv" + suffix + "Sofia");
    }
    return(x);
  }
  $(".template-authors").html(country_flag($(".template-authors").html()));
  $(".template-citation-authors").html(country_flag($(".template-citation-authors").html()));

  // add doctoral titles to authors
  function doct_tit(x) {
    if (typeof(x) != "undefined") {
      x = x.replace(/Author, maintainer/g, "Principal developer");
      x = x.replace(/Author, contributor/g, "Package maintainer");
      x = x.replace(/Thesis advisor/g, "(former) Doctoral advisor");
      // contributors
      x = x.replace("Alex", "Prof. Alex");
      x = x.replace("Andrew", "Dr. Andrew");
      x = x.replace("Annick", "Dr. Annick");
      x = x.replace("Anthony", "Dr. Anthony");
      x = x.replace("Bart", "Dr. Bart");
      x = x.replace("Bhanu", "Prof. Bhanu");
      x = x.replace("Casper", "Prof. Casper");
      x = x.replace("Christian", "Dr. Christian");
      x = x.replace("Corinna", "Dr. Corinna");
      x = x.replace("Dennis", "Dr. Dennis");
      x = x.replace("Gwen", "Prof. Gwen");
      x = x.replace("Jason", "Dr. Jason");
      x = x.replace("Javier", "Prof. Javier");
      x = x.replace("Jonas", "Dr. Jonas");
      x = x.replace("Judith", "Dr. Judith");
      x = x.replace("Larisse", "Dr. Larisse");
      x = x.replace("Matthijs", "Dr. Matthijs");
      x = x.replace("Peter", "Dr. Peter");
      x = x.replace("Rogier", "Dr. Rogier");
      x = x.replace("Sofia", "Dr. Sofia");
    }
    return(x);
  }
  $(".template-authors").html(doct_tit($(".template-authors").html()));
  $(".template-citation-authors").html(doct_tit($(".template-citation-authors").html()));
  $(".developers").html(doct_tit($(".developers").html()));
  $(".developers a[href='authors.html']").text("All contributors...");
});

$('head').append("<!-- Global site tag (gtag.js) - Google Analytics --> <script async src=\"https://www.googletagmanager.com/gtag/js?id=UA-172114740-1\"></script> <script> window.dataLayer = window.dataLayer || []; function gtag(){dataLayer.push(arguments);} gtag('js', new Date()); gtag('config', 'UA-172114740-1'); </script><!-- Matomo --><script type='text/javascript'> var _paq = _paq || []; /* tracker methods like 'setCustomDimension' should be called before 'trackPageView' */ _paq.push(['setDomains', ['*.msberends.github.io/AMR']]); _paq.push(['enableCrossDomainLinking']); _paq.push(['trackPageView']); _paq.push(['enableLinkTracking']); (function() { var u='https://analyse.uscloud.nl/'; _paq.push(['setTrackerUrl', u+'piwik.php']); _paq.push(['setSiteId', '3']); var d=document, g=d.createElement('script'), s=d.getElementsByTagName('script')[0]; g.type='text/javascript'; g.async=true; g.defer=true; g.src=u+'piwik.js'; s.parentNode.insertBefore(g,s);  })();</script><!-- End Matomo Code -->");
