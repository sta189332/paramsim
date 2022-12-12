library(pacman)
p_load(
  devtools,
  available,
  usethis,
  roxygen2,
  testthat,
  knitr,
  rmarkdown,
  rlang,
  sinew
)

devtools::has_devel()
pkgbuild::check_build_tools()
available::available("paramsim")

sinew::makeOxygen(arimasim)
usethis::use_git()
usethis::use_test()
devtools::test()
testthat::test_check(arimasimsearch)
usethis::use_cran_comments(open = rlang::is_interactive()) # before you submit to CRAN
devtools::release() # to send package to CRAN
