/*
# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #
*/

// Add updated Font Awesome 5.6.3 library
$('head').append('<!-- Updated Font Awesome library --><link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.6.3/css/all.css" integrity="sha384-UHRtZLI+pbxtHCWp1t77Bi1L4ZtiqrqD80Kn4Z8NTSRyMA2Fd33n5dQ8lWUE00s/" crossorigin="anonymous">');
// Add Disqus
$('head').append('<script src="https://amr-for-r.disqus.com/embed.js" data-timestamp="' + new Date()  + '"></script');

// Email template for new GitLab issues
//https://stackoverflow.com/a/33190494/4575331
//incoming+msberends-amr-9011429-5miwzuo1xo70wbz9r6fwv4dmg-issue@incoming.gitlab.com

// Edit footer
$( document ).ready(function() {

  // redirect to GitLab
  var url_old = window.location.href;
  var url_new = url_old.replace("github", "gitlab");
  if (url_old != url_new) {
    window.location.replace(url_new);
  }

  if ($(".template-article").length > 0) {
    // PR for 'R for Data Science' on How To pages
    $('#sidebar').prepend(
    '<div id="r4ds">' +
    '  <a target="_blank" href="https://r4ds.had.co.nz/">' +
    '    Learn R reading this great book: R for Data Science.' +
    '    <br><br>' +
    '    Click to read it online - it was published for free.' +
    '    <img src="https://gitlab.com/msberends/AMR/raw/master/docs/cover_r4ds.png" height="100px">' +
    '  </a>     ' +
    '  <hr>' +
    '</div>');
  }

  // Add Disqus to all pages
  var disqus =
    '<button class="btn btn-primary" type="button" data-toggle="collapse" data-target="#collapseDisqus" aria-expanded="false" aria-controls="collapseDisqus">' +
      '  Having a question? Or knowing something we don\'t? Click here to open comments.' +
    '</button>' +
    '<div class="collapse" id="collapseDisqus">' +
      '<div id="disqus_thread"></div>' +
    '</div>';

  $(disqus).insertBefore('footer');
  $('#disqus_thread footer').remove();

  // Alter footer
  $('footer').html(
    '<div>' +
      '<p>' + $('footer .copyright p').html().replace(
        "Developed by",
        '<code>AMR</code> (for R). Developed at the <a href="https://www.rug.nl">University of Groningen</a>.<br>Authors:') + '</p>' +
      '<a href="https://www.rug.nl"><img src="https://gitlab.com/msberends/AMR/raw/master/docs/logo_rug.png" class="footer_logo"></a>' +
    '</div>');
  // all links should open in new tab/window
  $('footer').html($('footer').html().replace(/href/g, 'target="_blank" href'));

  // doctoral titles of authors
  function doct_tit(x) {
    if (typeof(x) != "undefined") {
      x = x.replace("Corinna", "Dr Corinna");
      x = x.replace("Alex", "Prof Dr Alex");
      x = x.replace("Bhanu", "Prof Dr Bhanu");
      x = x.replace(/Author, thesis advisor/g, "Doctoral advisor");
      x = x.replace(/Authors/g, "aut_plural");
      x = x.replace(/Author, maintainer[.]?/g, "");
      x = x.replace(/Author[.]?/g, "");
      x = x.replace(/aut_plural/g, "Authors");
    }
    return(x);
  }
  $(".template-authors").html(doct_tit($(".template-authors").html()));
  $(".developers").html(doct_tit($(".developers").html()));
  $("footer").html(doct_tit($("footer").html()));

  // Edit title of manual
  $('.template-reference-index h1').text('Manual');
});

$('head').append("<!-- Matomo --><script type='text/javascript'> var _paq = _paq || []; /* tracker methods like 'setCustomDimension' should be called before 'trackPageView' */ _paq.push(['setDomains', ['*.msberends.gitlab.io/AMR']]); _paq.push(['enableCrossDomainLinking']); _paq.push(['trackPageView']); _paq.push(['enableLinkTracking']); (function() { var u='https://analyse.uscloud.nl/'; _paq.push(['setTrackerUrl', u+'piwik.php']); _paq.push(['setSiteId', '3']); var d=document, g=d.createElement('script'), s=d.getElementsByTagName('script')[0]; g.type='text/javascript'; g.async=true; g.defer=true; g.src=u+'piwik.js'; s.parentNode.insertBefore(g,s);  })();</script><!-- End Matomo Code -->");
