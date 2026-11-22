# Builds the pkgdown site of nmpilot, with the task guide on the home page.
# Run from the package root:   source("pkgdown/build_site.R")
# Needs the pkgdown and pkgload packages (pkgload comes with devtools).
#
# 1. writes the guide (nm_guide()) to pkgdown/assets/guide.html, copied by
#    pkgdown to the root of the site;
# 2. writes the home page pkgdown/index.md: README.md followed by the guide,
#    shown in a frame (README.md itself is not changed);
# 3. builds the site in docs/.
# pkgdown/index.md and pkgdown/assets/guide.html are generated: do not edit
# them, edit README.md and inst/guide/guide_tasks.csv instead.

pkgload::load_all(".", quiet = TRUE)   # nm_guide() reads inst/guide/guide_tasks.csv

dir.create("pkgdown/assets", recursive = TRUE, showWarnings = FALSE)
nm_guide(file = "pkgdown/assets/guide.html", open = FALSE)

readme <- readLines("README.md", encoding = "UTF-8")
home <- c(
  readme,
  "",
  "## Task guide",
  "",
  paste("Choose a stage of the analysis, then a task, to see the matching",
        "nmpilot function and an example of R code. The same guide opens in R",
        "with `nm_guide()` and in RStudio from the Addins menu",
        "([full page](guide.html))."),
  "",
  paste0('<iframe src="guide.html" title="nmpilot task guide" loading="lazy" ',
         'style="width:100%;height:80vh;border:1px solid #d8dee4;',
         'border-radius:8px;"></iframe>'),
  "")
writeLines(home, "pkgdown/index.md", useBytes = TRUE)

pkgdown::build_site(preview = interactive())
