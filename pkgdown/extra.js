/*
# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
#                                                                      #
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
*/

$(document).ready(function() {

  // Replace 'Developers' with 'Maintainers' on the main page
  $(".developers h2").text("Maintainers");

  // replace \donttest and \dontrun texts in Examples
  if ($("#ref-examples ~ div pre").length > 0) {
    $("#ref-examples ~ div pre").html($("#ref-examples ~ div pre").html().replaceAll("# \\donttest{", ""));
    $("#ref-examples ~ div pre").html($("#ref-examples ~ div pre").html().replaceAll("# \\dontrun{", ""));
    $("#ref-examples ~ div pre").html($("#ref-examples ~ div pre").html().replaceAll("# }", ""));
  }
  
  // remove leading newline in code examples on changelog
  if ($("body .template-news").length > 0) {
    $("body .template-news").html($("body .template-news").html().replaceAll('sourceCode R">\n<span', 'sourceCode R"><span'));
    $("body .template-news").html($("body .template-news").html().replaceAll('sourceCode R"><span></span>\n<span', 'sourceCode R"><span'));
  }

  // add doctoral titles to authors
  function doct_tit(x) {
    if (typeof(x) != "undefined") {
      x = x.replace(/Author, maintainer/g, "Principle developer");
      x = x.replace(/Author, contributor/g, "Contributing maintainer");
      x = x.replace(/Author, thesis advisor/g, "Doctoral advisor");
      x = x.replace(/Thesis advisor/g, "Doctoral advisor");
      // authors
      x = x.replace("Alex", "Prof. Dr. Alex");
      x = x.replace("Anthony", "Dr. Anthony");
      x = x.replace("Bart", "Dr. Bart");
      x = x.replace("Bhanu", "Prof. Dr. Bhanu");
      x = x.replace("Casper", "Prof. Dr. Casper");
      x = x.replace("Christian", "Dr. Christian");
      x = x.replace("Corinna", "Dr. Corinna");
      x = x.replace("Dennis", "Dr. Dennis");
      x = x.replace("Gwen", "Dr. Gwen");
      x = x.replace("Judith", "Dr. Judith");
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
