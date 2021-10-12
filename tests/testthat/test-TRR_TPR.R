context("Test TRR and TPR with .sim function, bat and square.")

testthat::skip('skip')
RNGkind("L'Ecuyer-CMRG")

test_that("TRR works with .sim function", {
  # Set digits option
  options(digits = 4)

  # # Install the package if it is not installed yet
  # install.packages("TRES")
  library("TRES")

  ## Set up the kind of random number generator (RNG)
  RNGkind("L'Ecuyer-CMRG")
  ## ----------------------------- Section 3.1 ----------------------------- ##

  # The TRR model: estimation, coefficient plots, estimation error of coefficient and subspace distance
  # Load dataset "bat"
  data("bat", package = "TRES")
  str(bat)

  # Fitting the TRR model with different methods
  fit_ols1 <- TRR.fit(bat$x, bat$y, method = "standard")
  fit_1d1 <- TRR.fit(bat$x, bat$y, u = c(14, 14), method = "1D")
  fit_pls1 <- TRR.fit(bat$x, bat$y, u = c(14, 14), method = "PLS")

  print(fit_1d1)
})
