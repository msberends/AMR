// Add updated Font Awesome 5.6.3 library
$('head').append('<!-- Updated Font Awesome library --><link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.6.3/css/all.css" integrity="sha384-UHRtZLI+pbxtHCWp1t77Bi1L4ZtiqrqD80Kn4Z8NTSRyMA2Fd33n5dQ8lWUE00s/" crossorigin="anonymous">');

/* edit footer */
$( document ).ready(function() {
  $('footer').html('<p>' +
                   $('footer .copyright p').html().replace("Developed by",
                                                           "<code>AMR</code> (for R). Developed at the University of Groningen.<br>Authors:") +
                   '</p>');
//$('footer').prepend("<div class='university'/>");
});
