testthat::context("Included Tests from Respec ")

ql <- c(68.0, 20.0, 345.0, 139.0, 680.0, 970.0, 355.0, 76.0, 110.0, 12.0, 
        117.0, 102.0, 74.0, 82.0, 221.0, 84.0, 55.0, 117.0, 225.0, 94.0, 
        0.0, 320.0, 0.0, 0.0, 41.0, 27.0, 205.0, 174.0, 100.0, 108.0, 
        387.0, 29.0, 111.0, 178.0)
qu <- c(68.0, 20.0, 345.0, 139.0, 680.0, 970.0, 355.0, 76.0, 110.0, 12.0, 
        117.0, 102.0, 74.0, 82.0, 221.0, 84.0, 55.0, 117.0, 225.0, 94.0, 
        40.0, 320.0, 40.0, 40.0, 41.0, 27.0, 205.0, 174.0, 100.0, 108.0, 
        387.0, 29.0, 111.0, 178.0)

ql_truth <- c(-0.69, 2.04, 0.36, 0.17, -0.16, 0.59, -0.84, -0.63, 1.03, -0.35, 
              0.86, 0.69, 0.81, -0.57, -0.47, 1.24, 0.29, -0.33, 0.22, 0.39, 
              0.08, -0.12, 1.13, -0.25, 0.67, -0.28, 0.17, -0.23, -0.59, -0.7, 
              -0.45, -0.26, -0.56, -0.47, -0.26, 0.26, 0.35, 0.64, 0.11, -0.54, 
              0, 0.22, 0.03, -0.12, -0.06, -0.08, 0.51, -0.58, 0.17, 0.01, 
              -0.76, -1.1, 0.25, 0.13, -0.63, -0.53, 1.14, 1.05, -0.42, 0.36, 
              0.52, -0.54, 0.27, 0.86, -0.3, -0.37, 0.55, 0.24, -0.05, -0.38, 
              -0.04, 0.55, -0.62, 0.26, -0.28, -0.15, -0.13, -0.93, 0.46, 0.52, 
              -0.34, -0.15, -0.08, 0.04, 0.42, 0.63, 0.66, -0.33, 0.05, -0.14, 
              -0.34, 0.76, 0.36, 0.93, 0.27, 0.08, 0.22, 0.93, -0.67, 0.3, 0.08,
              0.65, -0.28, 0.38, -0.05, 0.17, 0.65, -0.07, -0.29, -0.33, -0.29, 
              0.4, 0.61, -0.42, -0.14, -0.65, 0.41, 0.4, 0.19, -0.67, 0.2, -0.41, 
              -0.28, 1.12, 1.3, -0.61, 0.44, 1.07, 0.14, 0.45, -0.77, -0.18, 
              -0.03, -0.05, -0.31, -0.54, 0.45, 0.54, 0.8, -0.24, 0.73, 0.45, 
              -0.63, -0.33, 0.09, -0.25, 0.37, 0.87, 0.43, -0.86, 0.18, -0.3, 
              -0.56, -0.15, -0.73, 0.41, -0.47, 0.3, 0.34, 1.4, -0.95, -0.49, 
              1.27, 0.34, 0.2, -0.74, 0.36, 0.57, -0.15, -0.53, 0.73, 0.27, 
              0.06, 0.34, 1.46, -0.48, 1.26, 0.33, -0.08, -0.45, 0.52, -0.84, 
              0.55, -0.25, -0.36, 0.02, -0.15, -0.53, 0.46, 0, 0.38, -0.57,
              0.23, -0.59, 0.34, -0.76, 0.09, -0.3, 0.24, 0.92)

qu_truth <- c(-0.69, 2.04, 0.36, 0.17, -0.16, 0.59, -0.84, -0.63, 1.03, -0.35,
                 0.86, 0.69, 0.81, -0.57, -0.47, 1.24, 0.29, -0.33, 0.22, 0.39,
                 0.08, -0.12, 1.13, -0.25, 0.67, -0.28, 0.17, -0.23, -0.59, -0.7,
                 -0.45, -0.26, -0.56, -0.47, -0.26, 0.26, 0.35, 0.64, 0.11, -0.54,
                 0, 0.22, 0.03, -0.12, -0.06, -0.08, 0.51, -0.58, 0.17, 0.01, -0.76,
                 -1.1, 0.25, 0.13, -0.63, -0.53, 1.14, 1.05, -0.42, 0.36, 0.52,
                 -0.54, 0.27, 0.86, -0.3, -0.37, 0.55, 0.24, -0.05, -0.38, -0.04,
                 0.55, -0.62, 0.26, -0.28, -0.15, -0.13, -0.93, 0.46, 0.52, -0.34,
                 -0.15, -0.08, 0.04, 0.42, 0.63, 0.66, -0.33, 0.05, -0.14, -0.34,
                 0.76, 0.36, 0.93, 0.27, 0.08, 0.22, 0.93, -0.67, 0.3, 0.08, 0.65,
                 -0.28, 0.38, -0.05, 0.17, 0.65, -0.07, -0.29, -0.33, -0.29, 0.4,
                 0.61, -0.42, -0.14, -0.65, 0.41, 0.4, 0.19, -0.67, 0.2, -0.41,
                 -0.28, 1.12, 1.3, -0.61, 0.44, 1.07, 0.14, 0.45, -0.77, -0.18,
                 -0.03, -0.05, -0.31, -0.54, 0.45, 0.54, 0.8, -0.24, 0.73, 0.45,
                 -0.63, -0.33, 0.09, -0.25, 0.37, 0.87, 0.43, -0.86, 0.18, -0.3,
                 -0.56, -0.15, -0.73, 0.41, -0.47, 0.3, 0.34, 1.4, -0.95, -0.49,
                 1.27, 0.34, 0.2, -0.74, 0.36, 0.57, -0.15, -0.53, 0.73, 0.27,
                 0.06, 0.34, 1.46, -0.48, 1.26, 0.33, -0.08, -0.45, 0.52, -0.84,
                 0.55, -0.25, -0.36, 0.02, -0.15, -0.53, 0.46, 0, 0.38, -0.57,
                 0.23, -0.59, 0.34, -0.76, 0.09, -0.3, 0.24, 0.92)

test_that("test_moms_p3", {

  #expected values from tests
  moms_x_orig = c(165.5409898, 40060.40734, 1.423055596)
  moms_x_ERL = c(165.5409898, 40060.40734, 1.385585543)
  moms_x_HWN = c(165.5409898, 40060.40734, 1.449544227)

  #input values for tests
  n = 34
  nG_orig = n * 1.0 * 0.160 / 0.206
  nG_ERL = n * 1.056 * 0.160 / 0.206
  nG_HWN = n * 0.962 * 0.160 / 0.206
  rG = -0.145
  mc_old = c(0.0, 1.0, 0.0)
  
  #running tests
  moms_orig <- moms_p3_R(ql,qu,rG,nG_orig,mc_old)
  expect_equal(moms_orig, moms_x_orig)
  
  moms_ERL <- moms_p3_R(ql,qu,rG,nG_ERL,mc_old)
  expect_equal(moms_ERL, moms_x_ERL)
  
  moms_HWN <- moms_p3_R(ql,qu,rG,nG_HWN,mc_old)
  expect_equal(moms_HWN, moms_x_HWN)
})

test_that("test_p3est_ema", {

  #expected values from tests
  moms_x_orig = c(167.236599216544, 39541.5780590362, 1.46729880561388)
  moms_x_ERL = c(167.238920642166, 39540.8926394831, 1.42946503056195)
  moms_x_HWN = c(167.234871013990, 39542.0884205924, 1.49402326308862)
  
  #input values for tests
  r_G = -0.145
  as_M_mse = 0.006306109073542849
  as_S2_mse = 0.002770797581276051
  as_G_mse = 0.155132529907670
  r_M = 0.0
  r_S2 = 1.0
  r_M_mse = -99.0
  r_S2_mse = -3996.0
  r_G_mse = 0.206115990877151
  
  #running tests 
  Wd = 1
  moms_orig <- p3est_ema_R (ql, qu, as_M_mse, as_S2_mse, as_G_mse, r_M, 
                         r_S2, r_G, r_M_mse, r_S2_mse, r_G_mse, Wd)
  expect_equal(moms_orig, moms_x_orig)
  
  Wd = 1.056
  moms_ERL <- p3est_ema_R (ql, qu, as_M_mse, as_S2_mse, as_G_mse, r_M, 
                            r_S2, r_G, r_M_mse, r_S2_mse, r_G_mse, Wd)
  expect_equal(moms_ERL, moms_x_ERL)
  
  Wd = 0.962
  moms_HWN <- p3est_ema_R (ql, qu, as_M_mse, as_S2_mse, as_G_mse, r_M, 
                            r_S2, r_G, r_M_mse, r_S2_mse, r_G_mse, Wd)
  expect_equal(moms_HWN, moms_x_HWN)
})

test_that("truth_test_moms_p3", {

  #expected values from tests
  moms_x = c(0.07415,	0.307708822,	0.452403784)

  #input values for tests
  nG = 0.0      
  rG = 0.0
  
  #running tests
  mc_old = c(0.0, 1.0, 0.0)
  moms = moms_p3_R(ql_truth,qu_truth,rG,nG,mc_old)
  expect_equal(moms, moms_x)
  
  #advanced test, with 50+ additional values
  #input values for tests
  ql_truth_2 <- append(ql_truth,rep(c(-10),50))
  qu_truth_2 <- append(qu_truth,rep(c(10),50))
  mc_old= c(moms[1], moms[2], moms[3])
  
  #running tests
  moms = moms_p3_R(ql_truth_2,qu_truth_2,rG,nG,mc_old)
  expect_equal(moms, moms_x)
  
})

test_that("truth_test_p3est_ema", {

    #expected values from tests
    moms_x = c(0.07415, 0.307708822,	0.452403784)

    #input values for tests
    n = 200
    Wd = 1.0
    r_G = 0.0         
    as_M_mse = 1    
    as_S2_mse = 1
    r_M = 0.0          
    r_S2 = 1.0         
    r_M_mse = -99
    r_S2_mse = -99
    r_G_mse = 1
    
    nthresh = 1
    nobs = n
    tl = -20
    tu = 20
    as_G_mse = mseg_all_R(nthresh, nobs, tl, tu, moms_x) 
    
    #running test 
    cmoms <- p3est_ema_R (ql_truth, qu_truth, as_M_mse, as_S2_mse, as_G_mse, r_M, 
                              r_S2, r_G, r_M_mse, r_S2_mse, r_G_mse, Wd)
    moms_x_new = (moms_x[3]*r_G_mse + r_G*as_G_mse)/(r_G_mse + as_G_mse)
    expect_equal(cmoms[1], moms_x[1])
    expect_equal(cmoms[2], moms_x[2])
    expect_equal(cmoms[3], moms_x_new)
    
    #advanced test, with 50+ additional values
    #input values for tests
    ql_truth_2 <- append(ql_truth,rep(c(-10),50))
    qu_truth_2 <- append(qu_truth,rep(c(10),50))
    n = 250
    nthresh = 2
    tl = c(-20, 10) 
    tu = c(20, 20)
    nobs = c(200 ,50)
    Wd = detrat(r_M, r_S2, moms_x, n, nthresh, nobs, tl, tu)
    r_G = 0.0          
    as_M_mse = 1   
    as_S2_mse = 1   
    as_G_mse = 0   
    r_M = 0.0         
    r_S2 = 1.0        
    r_M_mse = -99  
    r_S2_mse = -99
    r_G_mse = 1
    as_G_mse = mseg_all_R(nthresh,nobs,tl,tu,moms_x)   
    
    #running test
    cmoms <- p3est_ema_R (ql_truth_2, qu_truth_2, as_M_mse, as_S2_mse, as_G_mse, r_M, 
                          r_S2, r_G, r_M_mse, r_S2_mse, r_G_mse, Wd)
    moms_x_new = (moms_x[3]*r_G_mse + r_G*as_G_mse)/(r_G_mse + as_G_mse)
    expect_equal(cmoms[1], moms_x[1])
    expect_equal(cmoms[2], moms_x[2])
    expect_equal(cmoms[3], moms_x_new)
})

test_that("truth_test_var_mom", {
  
  thr = 1.0e-2
  
  #input values for tests
  nthresh = 1
  n_in = 200
  tl_in = -20.0
  tu_in = 20.0
  mc_in <- c(0.07415, 0.307708822, 0.452403784)

  #expected values from tests
  n = 200
  var_m_x = c(mc_in[2]/n,  2*mc_in[2]**2/(n-1), (6/n)*(1 + (9.0/6.0*(mc_in[3]**2)) + (15.0/48.0*(mc_in[3]**4))))
  
  #running tests 
  var_m <- var_mom_R(nthresh,n_in,tl_in,tu_in,mc_in)
  expect_equal(var_m[1,1], var_m_x[1], thr)
  expect_equal(var_m[2,2], var_m_x[2], thr)
  # expect_equal(var_m[3,3], var_m_x[3], thr) commented out in respec test 
  
  #advanced test, with 50+ additional values
  #input values for tests
  nthresh = 2
  tl_in = c(-20.0, 10)
  tu_in = c(20.0, 20)
  n_in = c(200, 50)
  #running tests 
  var_m <- var_mom_R(nthresh,n_in,tl_in,tu_in,mc_in)
  expect_equal(var_m[1,1], var_m_x[1], thr)
  expect_equal(var_m[2,2], var_m_x[2], thr) 
  # expect_equal(var_m[3,3], var_m_x[3], thr) commented out in respec test
  
})

test_that("truth_test_qP3", {
  
  # expected values from tests
  q_in = c(0.9980, 0.9900, 0.9000, 0.5000, 0.1000, 0.0100, 0.0020)
  
  #input values for tests
  p = c()
  g = 2.0
  mu <- 0.07415
  s2 <- 0.307708822
  alpha = 4/g**2
  beta = g/abs(g) * (s2/alpha)**0.5
  tau = mu - alpha * beta
  mc_in = c(mu, s2,g)
  
  # skew of +2
  for (i in 1:7) {
    p[i] = 1 - exp((tau - qP3_R(q_in[i],mc_in))/beta)
  }
  expect_equal(p, q_in)
  
  #input values for tests
  p = c()
  g = -2.0
  mu <- 0.07415
  s2 <- 0.307708822
  alpha = 4/g**2
  beta = g/abs(g) * (s2/alpha)**0.5
  tau = mu - alpha * beta
  mc_in = c(mu, s2,g)
  
  # repeat with skew of -2
  for (i in 1:7) {
    p[i] = exp((tau - qP3_R(q_in[i],mc_in))/beta)
  }
  expect_equal(p, q_in)
  
})

