/*
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
*/

// Add updated Font Awesome 5.8.2 library
$('head').append('<!-- Updated Font Awesome library --><link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.2/css/all.css">');

$(document).ready(function() {

  // remove version label from header
  $(".version.label").remove();

  // redirect GitLab to GitHub
  var url_old = window.location.href;
  var url_new = url_old.replace("gitlab", "github");
  if (url_old != url_new) {
    window.location.replace(url_new);
  }
  
  // Edit title of manual
  $('.template-reference-index h1').text('Manual');

  // replace 'Value' in manual with 'Returned value'
  $(".template-reference-topic h2#value").text("Returned value");
  
  // replace \donttest and \dontrun texts in Examples
  if ($(".ref-examples pre").length > 0) {
    $(".ref-examples pre").html($(".ref-examples pre").html().replaceAll("# \\donttest{", ""));
    $(".ref-examples pre").html($(".ref-examples pre").html().replaceAll("# \\dontrun{", ""));
    $(".ref-examples pre").html($(".ref-examples pre").html().replaceAll("# }", ""));
  }
  
  // remove leading newline in code examples on changelog
  if ($("body .template-news").length > 0) {
    $("body .template-news").html($("body .template-news").html().replaceAll('sourceCode R">\n<span', 'sourceCode R"><span'));
  }
  // change H1  header on dev version on changelog, since pkgdown uses the version number from the installed version
  // (rather then using the DESCRIPTION file)
  $("h1[id^=unreleased]").text("Current development version");
  
  // PR for 'R for Data Science' on How To pages
  if ($(".template-article").length > 0) {
    $('#pkgdown-sidebar').prepend(
    '<div id="r4ds">' +
    '  <div class="txt">' +
    '    <p>' +
    '      Learn R reading this great book: <i>R for Data Science</i>.' +
    '      <br><br>' +
    '      <a target="_blank" href="https://r4ds.had.co.nz/">' +
    '        Click to read it online - it was published for free.' +
    '      </a>' +
    '    </p>' +
    '  </div>' +
    '  <div class="img">' +
    '    <a target="_blank" href="https://r4ds.had.co.nz/">' +
    '      <img src="https://github.com/msberends/AMR/raw/main/docs/cover_r4ds.png">' +
    '    </a>' +
    '  </div>' +
    '</div>');
  }

  // edit footer
  $('footer').html(
    '<div>' +
      '<p><code>AMR</code> (for R). Developed at the <a target="_blank" href="https://www.rug.nl">University of Groningen</a> in collaboration with non-profit organisations<br><a target="_blank" href="https://www.certe.nl">Certe Medical Diagnostics and Advice Foundation</a> and <a target="_blank" href="https://www.umcg.nl">University Medical Center Groningen</a>.</p>' +
            '<a target="_blank" href="https://www.rug.nl"><img src="https://github.com/msberends/AMR/raw/main/docs/logo_rug.png" class="footer_logo"></a>' + 
    '</div>');

  // doctoral titles of authors
  function doct_tit(x) {
    if (typeof(x) != "undefined") {
      // authors
      x = x.replace(/Author, maintainer/g, "Maintainer");
      x = x.replace(/Author, contributor/g, "Maintainer");
      x = x.replace(/Author, thesis advisor/g, "Doctoral advisor");
      x = x.replace(/Thesis advisor/g, "Doctoral advisor");
      x = x.replace("Matthijs", "Dr Matthijs");
      x = x.replace("Christian", "Dr Christian");
      x = x.replace("Alex", "Prof. Alex");
      x = x.replace("Bhanu", "Prof. Bhanu");
      x = x.replace("Casper", "Prof. Casper");
      x = x.replace("Corinna", "Dr Corinna");
      // others
      x = x.replace("Bart", "Dr Bart");
      x = x.replace("Sofia", "Dr Sofia");
      x = x.replace("Dennis", "Dr Dennis");
      x = x.replace("Judith", "Dr Judith");
      x = x.replace("Gwen", "Dr Gwen");
      x = x.replace("Anthony", "Dr Anthony");
      x = x.replace("Rogier", "Dr Rogier");
    }
    return(x);
  }
  $(".template-authors").html(doct_tit($(".template-authors").html()));
  $(".template-citation-authors").html(doct_tit($(".template-citation-authors").html()));
  $('.template-citation-authors h1').eq(0).html("All contributors of the <code>AMR</code> package");
  $('.template-citation-authors h1').eq(1).html("How to cite the <code>AMR</code> package");
  $(".developers").html(doct_tit($(".developers").html()));
  $(".developers a[href='authors.html']").text("All contributors...");
});

$('head').append("<!-- Global site tag (gtag.js) - Google Analytics --> <script async src=\"https://www.googletagmanager.com/gtag/js?id=UA-172114740-1\"></script> <script> window.dataLayer = window.dataLayer || []; function gtag(){dataLayer.push(arguments);} gtag('js', new Date()); gtag('config', 'UA-172114740-1'); </script><!-- Matomo --><script type='text/javascript'> var _paq = _paq || []; /* tracker methods like 'setCustomDimension' should be called before 'trackPageView' */ _paq.push(['setDomains', ['*.msberends.github.io/AMR']]); _paq.push(['enableCrossDomainLinking']); _paq.push(['trackPageView']); _paq.push(['enableLinkTracking']); (function() { var u='https://analyse.uscloud.nl/'; _paq.push(['setTrackerUrl', u+'piwik.php']); _paq.push(['setSiteId', '3']); var d=document, g=d.createElement('script'), s=d.getElementsByTagName('script')[0]; g.type='text/javascript'; g.async=true; g.defer=true; g.src=u+'piwik.js'; s.parentNode.insertBefore(g,s);  })();</script><!-- End Matomo Code -->");
