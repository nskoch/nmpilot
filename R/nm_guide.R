#' Interactive task guide
#'
#' Builds and opens an HTML page that lists the steps of a population PK
#' analysis with nmpilot. Buttons group the tasks by stage (data preparation,
#' exploration, control streams, running, results, diagnostics, model
#' selection, validation, reporting); selecting a task shows the matching
#' nmpilot function(s), a short description and an example of R code that can
#' be copied. Tasks whose function is not available yet are shown as
#' "planned".
#'
#' The page is generated from the file `guide/guide_tasks.csv` installed with
#' the package (one row per task: stage, task, functions, description, code,
#' status), so that it follows the package. It is a single self-contained
#' file: no internet connection is needed. The examples use the example
#' project of the package, copied to a temporary folder (see the "Setup" block
#' of the page); `nmfe` must be set to the path of your NONMEM launcher.
#'
#' @param output_file Output path of the HTML page. Default: `nmpilot_guide.html` in
#'   the temporary directory of the R session. (Formerly `file`, still accepted with a warning.)
#' @param ... Former argument name, accepted with a warning until version 1.0.
#' @param open Open the page in the web browser (default: in interactive
#'   sessions only).
#'
#' @return The path of the HTML file, invisibly.
#' @export
#'
#' @examples
#' p <- nm_guide(open = FALSE)
#' file.exists(p)
#' \dontrun{
#' nm_guide()     # opens the guide in the browser
#' }
nm_guide <- function(output_file = NULL, open = interactive(), ...) {
  dep <- .nm_deprecated_args(list(...), c(file = "output_file"), "nm_guide")
  if (!is.null(dep$output_file)) output_file <- dep$output_file
  file <- output_file
  tasks <- .nm_guide_tasks()
  version <- tryCatch(as.character(utils::packageVersion("nmpilot")),
                      error = function(e) "")
  html <- .nm_guide_html(tasks, version = version)
  if (is.null(file)) file <- file.path(tempdir(), "nmpilot_guide.html")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  con <- base::file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(html, con)
  if (isTRUE(open)) utils::browseURL(normalizePath(file, winslash = "/"))
  invisible(file)
}

#' RStudio addin: open the task guide
#'
#' Called by the RStudio addin "nmpilot task guide" (menu Addins). Writes the
#' guide with [nm_guide()] in the temporary directory and opens it in the
#' RStudio Viewer pane, or in the web browser outside RStudio.
#'
#' @return The path of the HTML file, invisibly.
#' @keywords internal
#' @export
nm_guide_addin <- function() {
  file <- nm_guide(open = FALSE)
  viewer <- getOption("viewer")
  if (is.function(viewer)) {
    viewer(file)
  } else {
    utils::browseURL(normalizePath(file, winslash = "/"))
  }
  invisible(file)
}

# Reads and checks the task table of the guide.
.nm_guide_tasks <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("guide", "guide_tasks.csv", package = "nmpilot")
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("The guide data file (guide/guide_tasks.csv) was not found. ",
         "Reinstall nmpilot.", call. = FALSE)
  }
  tasks <- utils::read.csv(path, stringsAsFactors = FALSE,
                           fileEncoding = "UTF-8", na.strings = character(0))
  needed <- c("task_id", "stage_id", "stage", "task", "functions",
              "description", "code", "status")
  missing <- setdiff(needed, names(tasks))
  if (length(missing)) {
    stop("Columns missing from guide_tasks.csv: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(unique(tasks$status), c("available", "planned"))
  if (length(bad)) {
    stop("Unknown status in guide_tasks.csv: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  tasks
}

.nm_guide_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

# Builds the HTML page (character vector of lines) from the task table.
.nm_guide_html <- function(tasks, version = "") {
  esc <- .nm_guide_escape
  stage_ids <- unique(tasks$stage_id)
  stage_lab <- tasks$stage[match(stage_ids, tasks$stage_id)]

  setup <- paste(
    "library(nmpilot)",
    "# Copy the example project to a temporary folder (runs write files there)",
    "file.copy(system.file(\"extdata\", \"example\", package = \"nmpilot\"),",
    "          tempdir(), recursive = TRUE)",
    "ex   <- file.path(tempdir(), \"example\")",
    "nmfe <- \"C:/nm76g64/run/nmfe76.bat\"   # path to your nmfe launcher",
    sep = "\n")

  stage_buttons <- sprintf(
    '<button class="stage" data-stage="%s">%s <span class="n">%d</span></button>',
    esc(stage_ids), esc(stage_lab),
    vapply(stage_ids, function(s) sum(tasks$stage_id == s), integer(1)))

  task_buttons <- sprintf(
    '<button class="task" data-stage="%s" data-task="%s">%s%s</button>',
    esc(tasks$stage_id), tasks$task_id, esc(tasks$task),
    ifelse(tasks$status == "planned", ' <span class="badge">planned</span>', ""))

  fun_chips <- function(f) {
    f <- trimws(strsplit(f, ",", fixed = TRUE)[[1]])
    f <- f[nzchar(f)]
    if (!length(f)) return("")
    paste(sprintf("<code class=\"fun\">%s</code>", esc(f)), collapse = " ")
  }
  help_line <- function(f) {
    f <- trimws(strsplit(f, ",", fixed = TRUE)[[1]])
    f <- f[nzchar(f)]
    if (!length(f)) return("")
    sprintf("<p class=\"help\">Help in R: %s</p>",
            paste(sprintf("<code>?%s</code>", esc(f)), collapse = ", "))
  }

  details <- vapply(seq_len(nrow(tasks)), function(i) {
    t <- tasks[i, ]
    planned <- t$status == "planned"
    body <- if (planned) {
      paste0("<p class=\"planned\">Not available yet: this task is planned ",
             "for a future version of nmpilot.</p>")
    } else {
      paste0(help_line(t$functions),
             "<div class=\"codebox\"><button class=\"copy\">Copy</button>",
             "<pre><code>", esc(t$code), "</code></pre></div>")
    }
    paste0(
      "<article class=\"detail\" data-task=\"", t$task_id, "\" hidden>",
      "<h2>", esc(t$task),
      if (planned) " <span class=\"badge\">planned</span>" else "", "</h2>",
      "<p class=\"funs\">", fun_chips(t$functions), "</p>",
      "<p class=\"desc\">", esc(t$description), "</p>",
      body, "</article>")
  }, character(1))

  css <- "
:root{--bg:#ffffff;--fg:#1f2328;--muted:#5d6670;--line:#d8dee4;--panel:#f6f8fa;
--accent:#1f6feb;--accent-fg:#ffffff;--code:#f3f4f6;--badge:#9a6700;--badge-bg:#fff8c5}
@media (prefers-color-scheme: dark){:root{--bg:#0d1117;--fg:#e6edf3;--muted:#8d96a0;
--line:#30363d;--panel:#161b22;--accent:#2f81f7;--accent-fg:#ffffff;--code:#1c2128;
--badge:#e3b341;--badge-bg:#3b2e00}}
*{box-sizing:border-box}
body{margin:0;font-family:system-ui,-apple-system,'Segoe UI',Roboto,sans-serif;
background:var(--bg);color:var(--fg);line-height:1.5}
header{padding:20px 28px 8px}
header h1{margin:0;font-size:1.5rem}
header p{margin:4px 0 0;color:var(--muted)}
main{padding:0 28px 28px}
details.setup{margin:12px 0 16px;border:1px solid var(--line);border-radius:8px;
background:var(--panel);padding:8px 12px}
details.setup summary{cursor:pointer;font-weight:600}
nav.stages{display:flex;flex-wrap:wrap;gap:8px;margin:8px 0 16px}
button{font:inherit;color:inherit}
button.stage{border:1px solid var(--line);background:var(--panel);border-radius:999px;
padding:6px 14px;cursor:pointer}
button.stage .n{color:var(--muted);font-size:.85em}
button.stage.active{background:var(--accent);border-color:var(--accent);color:var(--accent-fg)}
button.stage.active .n{color:var(--accent-fg)}
.layout{display:grid;grid-template-columns:minmax(220px,320px) 1fr;gap:20px}
@media (max-width:760px){.layout{grid-template-columns:1fr}}
.tasks{display:flex;flex-direction:column;gap:6px}
button.task{text-align:left;border:1px solid var(--line);background:var(--bg);
border-radius:6px;padding:8px 10px;cursor:pointer}
button.task:hover{border-color:var(--accent)}
button.task.active{border-color:var(--accent);box-shadow:inset 3px 0 0 var(--accent)}
.detail h2{margin:0 0 6px;font-size:1.2rem}
.desc{margin:8px 0}
.help{color:var(--muted);margin:8px 0}
code{font-family:ui-monospace,SFMono-Regular,Consolas,'Liberation Mono',monospace;font-size:.92em}
code.fun{background:var(--code);border:1px solid var(--line);border-radius:4px;padding:1px 6px}
.codebox{position:relative}
pre{background:var(--code);border:1px solid var(--line);border-radius:8px;padding:12px 14px;
overflow-x:auto;margin:8px 0}
button.copy{position:absolute;top:8px;right:8px;border:1px solid var(--line);
background:var(--panel);border-radius:6px;padding:2px 10px;cursor:pointer;font-size:.85em}
.badge{display:inline-block;font-size:.72rem;font-weight:600;color:var(--badge);
background:var(--badge-bg);border-radius:999px;padding:0 8px;vertical-align:middle}
.planned{color:var(--muted);font-style:italic}
footer{color:var(--muted);font-size:.85em;padding:0 28px 20px}
"

  js <- "
(function(){
  var stages = document.querySelectorAll('button.stage');
  var tasks = document.querySelectorAll('button.task');
  var details = document.querySelectorAll('article.detail');
  function showTask(id){
    tasks.forEach(function(b){ b.classList.toggle('active', b.dataset.task === id); });
    details.forEach(function(d){ d.hidden = d.dataset.task !== id; });
  }
  function showStage(s){
    stages.forEach(function(b){ b.classList.toggle('active', b.dataset.stage === s); });
    var first = null;
    tasks.forEach(function(b){
      var on = b.dataset.stage === s;
      b.hidden = !on;
      if (on && first === null) first = b.dataset.task;
    });
    showTask(first);
  }
  stages.forEach(function(b){ b.addEventListener('click', function(){ showStage(b.dataset.stage); }); });
  tasks.forEach(function(b){ b.addEventListener('click', function(){ showTask(b.dataset.task); }); });
  function fallbackCopy(text){
    var ta = document.createElement('textarea');
    ta.value = text; document.body.appendChild(ta); ta.select();
    try { document.execCommand('copy'); } catch (e) {}
    document.body.removeChild(ta);
  }
  document.querySelectorAll('button.copy').forEach(function(b){
    b.addEventListener('click', function(){
      var text = b.parentNode.querySelector('code').textContent;
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).catch(function(){ fallbackCopy(text); });
      } else { fallbackCopy(text); }
      b.textContent = 'Copied';
      setTimeout(function(){ b.textContent = 'Copy'; }, 1500);
    });
  });
  if (stages.length) showStage(stages[0].dataset.stage);
})();
"

  title <- if (nzchar(version)) paste0("nmpilot ", version, " \u2014 task guide")
           else "nmpilot \u2014 task guide"
  c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    paste0("<title>", esc(title), "</title>"),
    "<style>", css, "</style>",
    "</head>",
    "<body>",
    "<header>",
    paste0("<h1>", esc(title), "</h1>"),
    "<p>Choose a stage of the analysis, then a task, to see the nmpilot function and an example of R code.</p>",
    "</header>",
    "<main>",
    "<details class=\"setup\"><summary>Setup used by the examples</summary>",
    "<div class=\"codebox\"><button class=\"copy\">Copy</button>",
    paste0("<pre><code>", esc(setup), "</code></pre></div>"),
    "<p class=\"help\">Examples that only read results need <code>ex</code>; examples that run NONMEM also need <code>nmfe</code>. Objects such as <code>dat</code>, <code>s</code>, <code>runs</code> or <code>cmp</code> are created by earlier tasks of the guide.</p>",
    "</details>",
    "<nav class=\"stages\">", stage_buttons, "</nav>",
    "<div class=\"layout\">",
    "<div class=\"tasks\">", task_buttons, "</div>",
    "<div class=\"details\">", details, "</div>",
    "</div>",
    "</main>",
    paste0("<footer>Generated by nm_guide() from guide/guide_tasks.csv (",
           nrow(tasks), " tasks, ", sum(tasks$status == "planned"),
           " planned).</footer>"),
    "<script>", js, "</script>",
    "</body>",
    "</html>")
}
