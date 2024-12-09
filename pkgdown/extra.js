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
      x = x.replace("Alex", '<img src="lang_de.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Alex');
      x = x.replace("Andrew", '<img src="lang_us.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Andrew');
      x = x.replace("Anita", '<img src="lang_au.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Anita');
      x = x.replace("Annick", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Annick');
      x = x.replace("Anthony", '<img src="lang_en.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Anthony');
      x = x.replace("Anton", '<img src="lang_uk.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Anton');
      x = x.replace("Bart", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Bart');
      x = x.replace("Bhanu", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Bhanu');
      x = x.replace("Casper", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Casper');
      x = x.replace("Christian", '<img src="lang_de.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Christian');
      x = x.replace("Corinna", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Corinna');
      x = x.replace("Dennis", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Dennis');
      x = x.replace("Dmytro", '<img src="lang_uk.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Dmytro');
      x = x.replace("Eric", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Eric');
      x = x.replace("Erwin", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Erwin');
      x = x.replace("Gwen", '<img src="lang_en.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Gwen');
      x = x.replace("Jason", '<img src="lang_ca.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Jason');
      x = x.replace("Javier", '<img src="lang_ca.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Javier');
      x = x.replace("Jonas", '<img src="lang_de.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Jonas');
      x = x.replace("Judith", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Judith');
      x = x.replace("Larisse", '<img src="lang_sa.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Larisse');
      x = x.replace("Matthew", '<img src="lang_ca.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Matthew');
      x = x.replace("Matthijs", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Matthijs');
      x = x.replace("Peter", '<img src="lang_en.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Peter');
      x = x.replace("Rogier", '<img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Rogier');
      x = x.replace("Sofia", '<img src="lang_sv.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Sofia');
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
      x = x.replace("Gwen", "Dr. Gwen");
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
