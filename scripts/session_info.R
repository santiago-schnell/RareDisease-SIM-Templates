# session_info.R
# Captures R session information and writes it to session_info.txt.
# Run from the repository root to document the environment used for simulations.
#
# Usage:
#   source("scripts/session_info.R")
#   # or from the command line:
#   Rscript scripts/session_info.R

out_file <- "session_info.txt"

si <- sessionInfo()

# Write to file
sink(out_file)
cat("=== R Session Information ===\n")
cat("Captured:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(si)
sink()

cat("Session info written to:", out_file, "\n")
