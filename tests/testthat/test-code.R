context("Test reproducible examples in the paper.")

RNGkind("L'Ecuyer-CMRG")

test_that("Section 3.1", {
  data("bat")
  fit_ols1 <- TRR.fit(bat$x, bat$y, method = 'standard')
  fit_1d1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = "1D")
  fit_pls1 <- TRR.fit(bat$x, bat$y, u = c(14, 14), method = "PLS")
  dist_ols1 <- rTensor::fnorm(coef(fit_ols1) - bat$coefficients)
  dist_1d1 <- rTensor::fnorm(coef(fit_1d1) - bat$coefficients)
  dist_pls1 <- rTensor::fnorm(coef(fit_pls1) - bat$coefficients)
  Pdist_1d1 <- rep(NA_real_, 2)
  Pdist_pls1 <- rep(NA_real_, 2)
  for (i in 1:2) {
    Pdist_1d1[i] <- subspace(bat$Gamma[[i]], fit_1d1$Gamma[[i]])
    Pdist_pls1[i] <- subspace(bat$Gamma[[i]], fit_pls1$Gamma[[i]])
  }
  Pdist_1d1 <- sum(Pdist_1d1)
  Pdist_pls1 <- sum(Pdist_pls1)

  expect_equal(unname(c(summary(fit_ols1$coefficients@data))), c(-0.23776572,-0.008733723,-0.0011321,-0.0013459704,0.0028872783,0.12037262), tol = 1e-8)
  expect_equal(unname(c(summary(fit_1d1$coefficients@data))), c(-0.0090063837,-8.5386625e-05,-3.8275866e-07,0.0027569357,0.00012573236,0.092872601), tol = 1e-8)
  expect_equal(unname(c(summary(fit_pls1$coefficients@data))), c(-0.244712,-0.010421177,-0.0015230425,-0.0032793868,0.0022911425,0.13181997), tol = 1e-8)
  expect_equal(c(dist_ols1, dist_1d1, dist_pls1), c(0.972040026, 0.044169466, 1.191911596), tol = 1e-8)
  expect_equal(c(Pdist_1d1, Pdist_pls1), c(0.11357711, 1.46626919), tol = 1e-8)
})

test_that("Section 3.2", {
  data("square")
  fit_ols2 <- TPR.fit(square$x, square$y, method = "standard")
  fit_1d2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "1D")
  fit_pls2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "PLS")
  dist_ols2 <- rTensor::fnorm(coef(fit_ols2) - square$coefficients)
  dist_1d2 <- rTensor::fnorm(coef(fit_1d2) - square$coefficients)
  dist_pls2 <- rTensor::fnorm(coef(fit_pls2) - square$coefficients)

  Pdist_1d2 <- rep(NA_real_, 2)
  Pdist_pls2 <- rep(NA_real_, 2)
  for (i in 1:2) {
    Pdist_1d2[i] <- subspace(square$Gamma[[i]],fit_1d2$Gamma[[i]])
    Pdist_pls2[i] <- subspace(square$Gamma[[i]],fit_pls2$Gamma[[i]])
  }
  Pdist_1d2 <- sum(Pdist_1d2)
  Pdist_pls2 <- sum(Pdist_pls2)

  expect_equal(unname(c(summary(fit_ols2$coefficients@data))), c(-279.22534,-47.647935,0.3292953,0.1793506,50.11932,193.01954), tol = 1e-8)
  expect_equal(unname(c(summary(fit_1d2$coefficients@data))), c(0.0051724674,0.0064416012,0.041992609,0.072408052,0.06063525,0.48383079), tol = 1e-8)
  expect_equal(unname(c(summary(fit_pls2$coefficients@data))), c(0.0098859986,0.03624922,0.041738415,0.082750835,0.04797928,0.50897739), tol = 1e-8)
  expect_equal(c(dist_ols2, dist_1d2, dist_pls2), c(2225.6377647, 5.7978275,5.5905682), tol = 1e-8)
  expect_equal(c(Pdist_1d2, Pdist_pls2), c(1.41425441, 0.15631639), tol = 1e-8)
})


test_that("Section 3.3", {
  set.seed(1)
  # Dimension selection for both TRR and TPR models
  uhat1 <- TRRdim(bat$x, bat$y, maxdim = 32)
  uhat2 <- TPRdim(square$x, square$y, maxdim = 16)
  testthat::expect_equal(uhat1$bicval, c(-1226.3918, -1181.5243), tol = 1e-4)
  testthat::expect_equal(uhat1$mse, 1.9220775, tol = 1e-7)
  testthat::expect_equal(uhat2$mincv,  1.1054945, tol = 1e-7)
})

test_that("Section 3.5", {
  set.seed(1)
  data("EEG")
  # Dimension selection
  u_eeg <- TRRdim(EEG$x, EEG$y)
  testthat::expect_equal(u_eeg$bicval, c(1.8677789,1.6071015), tol = 1e-7)
})

test_that("Section 4.3", {
  # Compare the six core functions
  set.seed(1)
  p <- 20
  u <- 5

  ## Generate Gamma, M and U from Model (M1)
  data <- MenvU_sim(p, u, jitter = 1e-5)
  Gamma <- data$Gamma
  M <- data$M
  U <- data$U

  G <- vector("list", 8)
  G[[1]] <- simplsMU(M, U, u)
  G[[2]] <- ECD(M, U, u)
  G[[3]] <- manifold1D(M, U, u)
  G[[4]] <- OptM1D(M, U, u)
  G[[5]] <- manifoldFG(M, U, u)
  G[[6]] <- OptMFG(M, U, u)

  d <- rep(NA_real_, 8)
  for (i in 1:6) {
    d[i] <- subspace(G[[i]], Gamma)
  }
  # Compare the performance of the FG algorithm with different initial values: 1D estimator and randomly generated matrix.
  A <- matrix(runif(p*u), p, u)
  G[[7]] <- manifoldFG(M, U, u, Gamma_init = A)
  G[[8]] <- OptMFG(M, U, Gamma_init = A)
  for (i in 7:8) {
    d[i] <- subspace(G[[i]], Gamma)
  }
  expect_equal(d, c(1.2998376e-08,5.649722e-13,1.3689449e-07,6.13004e-13,1.315171e-07,6.5317157e-13,0.63245553,0.72155627), tol = 1e-8)
})


test_that("Section 4.4", {
  set.seed(1)
  p <- 50
  u <- 5
  n0 <- c(50, 70, 100, 200, 400, 800)
  uhat3 <- rep(NA_integer_, length(n0))
  for (i in seq_along(n0)) {
    n <- n0[i]
    data <- MenvU_sim(p, u, jitter = 1e-5, wishart = TRUE, n = n)
    M <- data$M
    U <- data$U
    output <- oneD_bic(M, U, n, maxdim = p/2)
    uhat3[i] <- output$u
  }
  expect_equal(uhat3, c(11,9,5,5,5,5))
})


