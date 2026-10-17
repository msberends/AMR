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
*/

$(function () {
  // add GPT assistant info
  $('aside').prepend('<div class="amr-gpt-assistant"><a target="_blank" href="https://chat.amr-for-r.org"><img src="https://amr-for-r.org/AMRforRGPT.svg"></a></div>');

  // split version number in navbar into main version and build number
  $('.nav-text').each(function () {
    const $el = $(this);
    const text = $.trim($el.text());
    const lastDotIndex = text.lastIndexOf('.');
    if (lastDotIndex > -1) {
      const main = text.substring(0, lastDotIndex);
      const build = text.substring(lastDotIndex);
      $el.html(`<span class="version-main">${main}</span><span class="version-build">${build}</span>`);
    }
  });

  // replace 'Developers' with 'Maintainers' on the main page, and "Contributors" on the Authors page
  $(".developers h2").text("Maintainers");
  $(".template-citation-authors h1:nth(0)").text("Contributors and Citation");
  $(".template-citation-authors h2:nth(0)").text("All contributors");
  $(".template-citation-authors h2:nth(1)").text("How to Cite");

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

  // reference index ('Manual'): collapse long function lists and filter topics, see pkgdown/templates/content-reference-index.html
  if ($(".amr-ref-item").length > 0) {
    const maxAliases = 8;
    const $body = $(".template-reference-index");
    const $input = $("#amr-ref-filter-input");
    const $status = $("#amr-ref-filter-status");
    const $items = $(".amr-ref-item");
    const pageText = {}; // full text per reference page, from search.json
    let pageTextLoaded = false;

    // keep the sticky filter field directly below the navbar, which hides when scrolling down
    const navbar = document.querySelector(".navbar.fixed-top");
    const filterBar = document.querySelector(".amr-ref-filter");
    function placeFilter() {
      const navBottom = navbar ? Math.max(0, navbar.getBoundingClientRect().bottom) : 0;
      $body.css("--amr-navbar-height", navBottom + "px");
      $(filterBar).toggleClass("amr-ref-stuck", filterBar.getBoundingClientRect().top <= navBottom + 1 && window.scrollY > 0);
    }
    let placing = false;
    $(window).on("scroll resize", function () {
      if (placing) return;
      placing = true;
      // follow the navbar during its slide transition
      const start = performance.now();
      (function step(now) {
        placeFilter();
        if (now - start < 400) requestAnimationFrame(step); else placing = false;
      })(start);
    });
    placeFilter();

    // headings that are jumped to (e.g. from 'On this page') must end up below the filter field, and below the navbar
    // when jumping up, since the navbar only shows when scrolling up
    $(document).on("click", 'a[href^="#"]', function () {
      const target = document.getElementById(decodeURIComponent($(this).attr("href").slice(1)));
      if (!target) return;
      const up = target.getBoundingClientRect().top < 0;
      document.documentElement.style.scrollPaddingTop = ((up && navbar ? navbar.offsetHeight : 0) + filterBar.offsetHeight + 8) + "px";
    });
    // the same when opening a link to a section, as the browser jumps before this script runs
    if (window.location.hash) {
      const target = document.getElementById(decodeURIComponent(window.location.hash.slice(1)));
      if (target) {
        document.documentElement.style.scrollPaddingTop = (filterBar.offsetHeight + 8) + "px";
        target.scrollIntoView();
      }
    }

    function updateMoreButton($item) {
      const nHidden = $item.find("code.amr-ref-extra:not(.amr-ref-match)").length;
      $item.find(".amr-ref-more").text("+ " + nHidden + " more").toggle(nHidden > 0);
    }

    $items.each(function () {
      const $item = $(this);
      const $codes = $item.find("dd code");
      // pkgdown puts backticks around non-syntactic names, e.g. `%like%`
      $codes.find("a").each(function () {
        const name = $(this).text();
        if (/^`.+`$/.test(name)) $(this).text(name.slice(1, -1));
      });
      $item.data("page", $item.find("dt a").attr("href"));
      $item.data("text", ($item.find("dt").text() + " " + $item.find("dd").text()).toLowerCase());
      // do not offer to show only one or two more
      if ($codes.length > maxAliases + 2) {
        $codes.slice(maxAliases).addClass("amr-ref-extra");
        $item.find("dd").append('<button type="button" class="amr-ref-more"></button>');
        updateMoreButton($item);
      }
    });
    let h2Title = "";
    $(".template-reference-index main").find("h2, h3").each(function () {
      if (this.tagName === "H2") h2Title = $(this).text().trim();
      else $(this).attr("data-amr-parent", h2Title);
    });
    $(".amr-ref-more").on("click", function () {
      $(this).closest(".amr-ref-item").addClass("amr-ref-expanded");
      $(this).remove();
    });

    function loadPageText() {
      if (pageTextLoaded) return;
      pageTextLoaded = true;
      $.getJSON($("#search-input").data("search-index") || "../search.json", function (data) {
        data.forEach(function (x) {
          if (x.dir !== "Reference") return;
          const page = x.path.split("/").pop();
          pageText[page] = (pageText[page] || "") + " " + (x.what + " " + x.text + " " + x.code).toLowerCase();
        });
        filter();
      });
    }

    function filter() {
      const terms = $.trim($input.val()).toLowerCase().split(/\s+/).filter(Boolean);
      $body.toggleClass("amr-ref-filtering", terms.length > 0);
      let nDirect = 0;
      let nText = 0;
      $items.each(function () {
        const $item = $(this);
        const $codes = $item.find("dd code");
        const direct = terms.length > 0 && terms.every(t => $item.data("text").includes(t));
        const text = terms.length > 0 && !direct && terms.every(t => (pageText[$item.data("page")] || "").includes(t));
        $item.toggleClass("amr-ref-hidden", terms.length > 0 && !direct && !text);
        $item.toggleClass("amr-ref-textmatch", text);
        // highlight the function names that match, which also reveals them if collapsed
        $codes.each(function () {
          const name = $(this).text().toLowerCase();
          $(this).toggleClass("amr-ref-match", direct && terms.some(t => name.includes(t)));
        });
        updateMoreButton($item);
        if (direct) nDirect++;
        if (text) nText++;
      });
      // hide headings (with their descriptions) that have no visible topic before the next heading of the same or a higher level;
      // this follows the document order, since pkgdown puts a heading and its topics in separate sibling blocks
      const nVisible = new Map();
      let h2 = null;
      let h3 = null;
      $(".template-reference-index main").find("h2, h3, .amr-ref-item").each(function () {
        if (this.tagName === "H2") {
          h2 = this;
          h3 = null;
          nVisible.set(this, 0);
        } else if (this.tagName === "H3") {
          h3 = this;
          nVisible.set(this, 0);
        } else if (!$(this).hasClass("amr-ref-hidden")) {
          if (h2) nVisible.set(h2, nVisible.get(h2) + 1);
          if (h3) nVisible.set(h3, nVisible.get(h3) + 1);
        }
      });
      nVisible.forEach((n, heading) => $(heading).closest(".section").toggleClass("amr-ref-hidden", n === 0));
      // blocks that only hold topics
      $(".template-reference-index main .section").not(":has(h2, h3)").each(function () {
        $(this).toggleClass("amr-ref-hidden", $(this).find(".amr-ref-item:not(.amr-ref-hidden)").length === 0);
      });
      const topics = n => n + (n === 1 ? " topic" : " topics");
      if (terms.length === 0) {
        $status.text("");
      } else if (nDirect + nText === 0) {
        $status.text(pageTextLoaded ? "No topics found" : "Searching...");
      } else if (nText === 0) {
        $status.text(topics(nDirect) + " found by name or title");
      } else if (nDirect === 0) {
        $status.text("No topics found by name or title, but " + topics(nText) + " mention it");
      } else {
        $status.text(topics(nDirect) + " found by name or title, and " + nText + " more mention it");
      }
      placeFilter();
    }

    $input.on("input", function () {
      loadPageText();
      filter();
    });
    $input.on("keydown", function (e) {
      if (e.key === "Escape") {
        $input.val("");
        filter();
      }
    });
    // press '/' anywhere to focus the filter field
    $(document).on("keydown", function (e) {
      if (e.key === "/" && !$(e.target).is("input, textarea")) {
        e.preventDefault();
        $input.trigger("focus");
      }
    });
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
    $('img[src="https://amr-for-r.org/AMRforRGPT.svg"]').attr('src', 'https://amr-for-r.org/AMRforRGPT_python.svg'); /* replace GPT logo */
  }
  
  // add country flags
  function country_flag(x) {
    if (typeof(x) != "undefined") {
      const prefix = '<img src="lang_';
      const suffix = '.svg" style="height: 13px !important; border: 1px solid #cccccc; aspect-ratio: 1 / 0.6; object-fit: cover; vertical-align: initial !important;"> ';
      // bit weird, but we refer to languages here! So "uk" is Ukrain, not UK - take "en" for that.
      x = x.replace("Aislinn", prefix + "en" + suffix + "Aislinn");
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
      x = x.replace("Jane", prefix + "au" + suffix + "Jane");
      x = x.replace("Jason", prefix + "ca" + suffix + "Jason");
      x = x.replace("Javier", prefix + "ca" + suffix + "Javier");
      x = x.replace("Jonas", prefix + "de" + suffix + "Jonas");
      x = x.replace("Judith", prefix + "nl" + suffix + "Judith");
      x = x.replace("Kathryn", prefix + "en" + suffix + "Kathryn");
      x = x.replace("Larisse", prefix + "sa" + suffix + "Larisse");
      x = x.replace("Matthew", prefix + "ca" + suffix + "Matthew");
      x = x.replace("Matthijs", prefix + "nl" + suffix + "Matthijs");
      x = x.replace("Natacha", prefix + "en" + suffix + "Natacha");
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
      x = x.replace(/Author, maintainer/g, "Project founder and lead developer");
      x = x.replace(/Author, contributor/g, "Package maintainer");
      x = x.replace(/Thesis advisor/g, "Former doctoral advisor");
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
      x = x.replace("Jane", "Dr. Jane");
      x = x.replace("Jason", "Dr. Jason");
      x = x.replace("Javier", "Prof. Javier");
      x = x.replace("Jonas", "Dr. Jonas");
      x = x.replace("Judith", "Dr. Judith");
      x = x.replace("Kathryn", "Prof. Kathryn");
      x = x.replace("Larisse", "Dr. Larisse");
      x = x.replace("Matthijs", "Dr. Matthijs");
      x = x.replace("Natacha", "Prof. Natacha");
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

  // add lines between collaborator types
  $('.template-citation-authors li:contains("Contributor")').first().before('<hr>');
  $('.template-citation-authors li:contains("Former doctoral advisor")').first().before('<hr>');

  // dark variants of figures: list a file name here and put a file with '_dark' before its extension next to it in pkgdown/assets,
  // e.g. "logo_umcg.svg" requires "logo_umcg_dark.svg"; images are swapped when the theme is dark, also after switching theme
  const darkVariants = [];
  function applyDarkVariants() {
    const dark = document.documentElement.getAttribute("data-bs-theme") === "dark";
    $("img").each(function () {
      const src = $(this).attr("src") || "";
      const file = src.split(/[?#]/)[0].split("/").pop();
      const base = file.replace(/_dark(\.\w+)$/, "$1");
      if (!darkVariants.includes(base)) return;
      const target = dark ? base.replace(/(\.\w+)$/, "_dark$1") : base;
      if (file !== target) $(this).attr("src", src.replace(file, target));
    });
  }
  applyDarkVariants();
  new MutationObserver(applyDarkVariants).observe(document.documentElement, { attributes: true, attributeFilter: ["data-bs-theme"] });
});

$('head').append("<!-- Matomo --><script type='text/javascript'> var _paq = _paq || []; /* tracker methods like 'setCustomDimension' should be called before 'trackPageView' */ _paq.push(['setDomains', ['*.amr-for-r.org']]); _paq.push(['enableCrossDomainLinking']); _paq.push(['trackPageView']); _paq.push(['enableLinkTracking']); (function() { var u='https://analyse.uscloud.nl/'; _paq.push(['setTrackerUrl', u+'piwik.php']); _paq.push(['setSiteId', '3']); var d=document, g=d.createElement('script'), s=d.getElementsByTagName('script')[0]; g.type='text/javascript'; g.async=true; g.defer=true; g.src=u+'piwik.js'; s.parentNode.insertBefore(g,s);  })();</script><!-- End Matomo Code -->");
