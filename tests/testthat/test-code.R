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

  expect_equal(unname(c(summary(fit_ols1$coefficients@data))), c(-0.237765722,-0.008733723,-0.001132100,-0.001345970,0.002887278,0.120372621), tol = 1e-8)
  expect_equal(unname(c(summary(fit_1d1$coefficients@data))), c(-8.968710e-03,-7.542956e-05,2.438858e-07,2.755897e-03,1.286767e-04,9.284800e-02), tol = 1e-8)
  expect_equal(unname(c(summary(fit_pls1$coefficients@data))), c(-0.2447121,-0.01042121,-0.00152305,-0.003279391,0.002291145,0.13182), tol = 1e-8)
  expect_equal(c(dist_ols1, dist_1d1, dist_pls1), c(0.97204003,0.044284581,1.1919118), tol = 1e-8)
  expect_equal(c(Pdist_1d1, Pdist_pls1), c(0.11434039, 1.46626921), tol = 1e-8)
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

  expect_equal(unname(c(summary(fit_ols2$coefficients@data))), c(-279.22535,-47.647966,0.32928287,0.17935058,50.119298,193.01952), tol = 1e-8)
  expect_equal(unname(c(summary(fit_1d2$coefficients@data))), c(0.0051790165,0.0064485306,0.042028368,0.072423178,0.060669237,0.48382142), tol = 1e-8)
  expect_equal(unname(c(summary(fit_pls1$coefficients@data))), c(0.0098860004,0.036249218,0.041738415,0.082750835,0.047979281,0.50897739), tol = 1e-8)
  expect_equal(c(dist_ols2, dist_1d2, dist_pls2), c(2225.6377950, 5.7978263, 5.5905682), tol = 1e-8)
  expect_equal(c(Pdist_1d2, Pdist_pls2), c(1.4142544,0.1563164), tol = 1e-8)
})


test_that("Section 3.3", {
  set.seed(1)
  # Dimension selection for both TRR and TPR models
  uhat1 <- TRRdim(bat$x, bat$y, maxdim = 32)
  uhat2 <- TPRdim(square$x, square$y, maxdim = 16)
  testthat::expect_equal(uhat1$bicval, c(-1261.6120, -1183.2608), tol = 1e-4)
  testthat::expect_equal(uhat1$mse, 1.9220851, tol = 1e-7)
  testthat::expect_equal(uhat2$mincv, 1.1054943, tol = 1e-7)
})

test_that("Section 3.5", {
  set.seed(1)
  data("EEG")
  # Dimension selection
  u_eeg <- TRRdim(EEG$x, EEG$y)
  testthat::expect_equal(u_eeg$bicval, c(1.8672561, 1.6075473), tol = 1e-7)
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
  expect_equal(d, c(1.4656091e-09,2.9683284e-13,5.1506432e-10,1.9246747e-09,4.5243985e-10,1.7448437e-09,0.4472136,0.60979972), tol = 1e-8)
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


