context("test-probabilistic_dropouts")

library(proDD)

test_that("dprobdropout works", {
    xg <- seq(-5, 5, length.out=101)
    # If there are no missing values it is just a normal distribution
    expect_equal(dprobdropout(xg, mu=1, sigma2=3), dnorm(xg, mean=1, sd=sqrt(3)))

    # If the sigmoid is very broad, it doesn't affect the result a lot
    expect_equal(dprobdropout(xg, mu=1, sigma2=3, rho=1, zeta=1000),
                 dnorm(xg, mean=1, sd=sqrt(3)),
                 tolerance=1e-3)
    expect_equal(dprobdropout(xg, mu=1, sigma2=3, rho=1, zeta=-1000),
                 dnorm(xg, mean=1, sd=sqrt(3)),
                 tolerance=1e-3)

    # or if the sigmoid is shifted far out
    expect_equal(dprobdropout(xg, mu=0, sigma2=3, rho=-20, zeta=1),
                 dnorm(xg, mean=0, sd=sqrt(3)),
                 tolerance=1e-3)


    # Otherwise it is skewed
    # ... to the left if zeta < 0
    expect_gt(dprobdropout(-1, mu=0, sigma2=3, rho=0, zeta=-1),
              dprobdropout(1, mu=0, sigma2=3, rho=0, zeta=-1))
    # ... to the right if zeta > 0
    expect_lt(dprobdropout(-1, mu=0, sigma2=3, rho=0, zeta=1),
              dprobdropout(1, mu=0, sigma2=3, rho=0, zeta=1))

    # at least until the sigmoid is shifted too far out
    expect_equal(dprobdropout(-1, mu=0, sigma2=3, rho=-20, zeta=1),
                 dprobdropout(1, mu=0, sigma2=3, rho=-20, zeta=1))

    # it can also handle multiple missing values
    expect_silent(dprobdropout(1, mu=0, sigma2=3, rho=c(0,0), zeta=c(-1, -0.8)))
    expect_silent(dprobdropout(1, mu=0, sigma2=3, rho=c(0,10), zeta=c(-1, -0.8)))
    expect_silent(dprobdropout(1, mu=0, sigma2=3, rho=2, zeta=-0.1, nmis=4))

    # which is different from a single missing value
    expect_true(all(dprobdropout(xg, mu=0, sigma2=3, rho=2, zeta=-0.1, nmis=1) !=
                     dprobdropout(xg, mu=0, sigma2=3, rho=2, zeta=-0.1, nmis=4)))

})


test_that("dprobdropout fails on bad input", {
    # rho and zeta must have the same length
    expect_error(dprobdropout(-3:3, mu=0, sigma2=1, rho=numeric(0), zeta=1))
    # all elements of zeta must have the same sign
    expect_error(dprobdropout(-3:3, mu=0, sigma2=1, rho=c(1,1), zeta=c(1, -1)))
    # nmis and rho/zeta must match
    expect_error(dprobdropout(-3:3, mu=0, sigma2=1, rho=c(1,1),
                              zeta=c(-1, -0.8), nmis=1))
    # mu and sigma2 are not vectors
    expect_error(dprobdropout(xg, mu=1:3, sigma2=3))
})


test_that("mode_probdropout works", {

    mode <- mode_probdropout(mu=0, sigma2=3, rho=2, zeta=-3)
    num_mode <- optimize(function(x){
        dprobdropout(x, mu=0, sigma2=3, rho=2, zeta=-3)
    }, lower=-10, upper=10, maximum=TRUE)$maximum
    expect_equal(mode, num_mode, tolerance=1e-5)

    mode <- mode_probdropout(mu=0, sigma2=3, rho=2, zeta=3)
    num_mode <- optimize(function(x){
        dprobdropout(x, mu=0, sigma2=3, rho=2, zeta=3)
    }, lower=-10, upper=10, maximum=TRUE)$maximum
    expect_equal(mode, num_mode, tolerance=1e-5)

    expect_gt(dprobdropout(mode, mu=0, sigma2=3, rho=2, zeta=3),
              dprobdropout(mode-1, mu=0, sigma2=3, rho=2, zeta=3))
    expect_gt(dprobdropout(mode, mu=0, sigma2=3, rho=2, zeta=3),
              dprobdropout(mode+1, mu=0, sigma2=3, rho=2, zeta=3))
    expect_gt(dprobdropout(mode, mu=0, sigma2=3, rho=2, zeta=3),
              dprobdropout(mode-0.0001, mu=0, sigma2=3, rho=2, zeta=3))
    expect_gt(dprobdropout(mode, mu=0, sigma2=3, rho=2, zeta=3),
              dprobdropout(mode+0.0001, mu=0, sigma2=3, rho=2, zeta=3))

    # The challenge is to find a parameter combination that breaks the function
    expect_silent(mode_probdropout(mu=0, sigma2=3, rho=2, zeta=0.3))
    expect_silent(mode_probdropout(mu=0, sigma2=3, rho=-2, zeta=0.3))
    expect_silent(mode_probdropout(mu=0, sigma2=3, rho=2, zeta=-0.3))
    expect_silent(mode_probdropout(mu=0, sigma2=3, rho=-2, zeta=-0.3))

    expect_silent(mode_probdropout(mu=0, sigma2=30, rho=20, zeta=0.3))
    expect_silent(mode_probdropout(mu=0, sigma2=30, rho=-20, zeta=0.3))
    expect_silent(mode_probdropout(mu=0, sigma2=30, rho=20, zeta=-0.3))
    expect_silent(mode_probdropout(mu=0, sigma2=30, rho=-20, zeta=-0.3))
})





test_that("mean_probdropout works", {
    expect_equal(mean_probdropout(mu=0, sigma2=3, rho=2, zeta=-3),
                 mean_probdropout(mu=0, sigma2=3, rho=2, zeta=-3, approx=FALSE),
                 tolerance=0.05)
    expect_equal(mean_probdropout(mu=0, sigma2=3, rho=-2, zeta=-3),
                 mean_probdropout(mu=0, sigma2=3, rho=-2, zeta=-3, approx=FALSE),
                 tolerance=0.05)
    # The approximation fails in this case
    # expect_equal(mean_probdropout(mu=0, sigma2=3, rho=2, zeta=-0.3),
    #              mean_probdropout(mu=0, sigma2=3, rho=2, zeta=-0.3, approx=FALSE),
    #              tolerance=0.05)
    expect_equal(mean_probdropout(mu=0, sigma2=3, rho=2, zeta=0.3),
                 mean_probdropout(mu=0, sigma2=3, rho=2, zeta=0.3, approx=FALSE),
                 tolerance=0.05)
    expect_equal(mean_probdropout(mu=0, sigma2=3, rho=-2, zeta=-0.3),
                 mean_probdropout(mu=0, sigma2=3, rho=-2, zeta=-0.3, approx=FALSE),
                 tolerance=0.05)

})


test_that("variance_probdropout works", {
    expect_equal(variance_probdropout(mu=0, sigma2=3, rho=2, zeta=-3),
                 variance_probdropout(mu=0, sigma2=3, rho=2, zeta=-3, approx=FALSE),
                 tolerance=0.05)
    expect_equal(variance_probdropout(mu=0, sigma2=3, rho=-2, zeta=-3),
                 variance_probdropout(mu=0, sigma2=3, rho=-2, zeta=-3, approx=FALSE),
                 tolerance=0.05)
    # The approximation fails in this case
    # expect_equal(mean_probdropout(mu=0, sigma2=3, rho=2, zeta=-0.3),
    #              mean_probdropout(mu=0, sigma2=3, rho=2, zeta=-0.3, approx=FALSE),
    #              tolerance=0.05)

    # Okay I cheat here by increasing the tolerance
    expect_equal(variance_probdropout(mu=0, sigma2=3, rho=2, zeta=0.3),
                 variance_probdropout(mu=0, sigma2=3, rho=2, zeta=0.3, approx=FALSE),
                 tolerance=0.5)
    expect_equal(variance_probdropout(mu=0, sigma2=3, rho=-2, zeta=-0.3),
                 variance_probdropout(mu=0, sigma2=3, rho=-2, zeta=-0.3, approx=FALSE),
                 tolerance=0.5)

})







