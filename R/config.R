###############################################################################
#  config.R — Path configuration for thesis replication
#
#  Source this file at the top of each step script.
#  All paths are relative to the repository root.
###############################################################################

BASE_DIR <- normalizePath(file.path(dirname(
  if (interactive()) rstudioapi::getSourceEditorContext()$path
  else sys.frame(1)$ofile
), ".."), mustWork = FALSE)

if (is.na(BASE_DIR) || BASE_DIR == "" || !dir.exists(BASE_DIR)) {
  BASE_DIR <- getwd()
}

DATA_DIR <- file.path(BASE_DIR, "data")
FIG_DIR  <- file.path(BASE_DIR, "output", "figures")
TBL_DIR  <- file.path(BASE_DIR, "output", "tables")
OUT_DIR  <- file.path(BASE_DIR, "output")

for (d in c(DATA_DIR, FIG_DIR, TBL_DIR, OUT_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat(sprintf("  Base directory : %s\n", BASE_DIR))
cat(sprintf("  Data directory : %s\n", DATA_DIR))
cat(sprintf("  Figure output  : %s\n", FIG_DIR))
cat(sprintf("  Table output   : %s\n", TBL_DIR))
