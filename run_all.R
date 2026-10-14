# Reproduce every result table (results/) and figure (figures/) in this repository.
#   Rscript run_all.R
# Run from the repository root. Total run time: about 2 minutes.
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
scripts <- sort(list.files("analysis", pattern = "^0[1-9].*\\.R$", full.names = TRUE))
for (s in scripts) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  source(s, local = new.env(), echo = FALSE)
  message(sprintf("     done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
message("\nAll analyses finished. Run the tests with: Rscript -e 'testthat::test_dir(\"tests/testthat\")'")
