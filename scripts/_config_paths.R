# scripts/_config_paths.R
# Single source of truth for portable paths.

get_root <- function() {
  # Priority: explicit env var -> --root arg in caller -> getwd()
  root <- Sys.getenv("PIPELINE_ROOT", unset = NA_character_)
  if (!is.na(root) && nzchar(root)) return(normalizePath(root, mustWork = TRUE))
  normalizePath(getwd(), mustWork = TRUE)
}

ROOT <- get_root()

DIRS <- list(
  root    = ROOT,
  data    = file.path(ROOT, "data"),
  gnomad  = file.path(ROOT, "data", "gnomad"),
  clinvar = file.path(ROOT, "data", "clinvar"),
  pext    = file.path(ROOT, "data", "pext"),
  domains = file.path(ROOT, "data", "domains"),
  tables  = file.path(ROOT, "tables"),
  figures = file.path(ROOT, "figures"),
  cache   = file.path(ROOT, "cache"),
  logs    = file.path(ROOT, "logs")
)

invisible(lapply(DIRS, dir.create, showWarnings = FALSE, recursive = TRUE))

p <- function(...) file.path(ROOT, ...)
pt <- function(...) file.path(DIRS$tables, ...)
pf <- function(...) file.path(DIRS$figures, ...)
pc <- function(...) file.path(DIRS$cache, ...)
pd <- function(...) file.path(DIRS$data, ...)
