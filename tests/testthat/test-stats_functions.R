context("test-stats_function")

library(proDD)

test_that("invprobit works", {
    xg <- seq(-10, 10, length.out = 1001)
    expect_equal(invprobit(xg, rho=0, zeta=1), pnorm(xg, mean=0, sd=1))
    expect_equal(invprobit(xg, rho=seq(-1, 1, length.out=1001), zeta=1),
                 pnorm(xg, mean=seq(-1, 1, length.out=1001), sd=1))
    expect_equal(invprobit(xg, rho=0, zeta=-1),
                 pnorm(xg, mean=0, sd=1, lower.tail = FALSE))

    # check log
    expect_equal(invprobit(xg, rho=0, zeta=-1, log=TRUE),
                 pnorm(xg, mean=0, sd=1, lower.tail = FALSE, log.p=TRUE))
    expect_equal(invprobit(xg, rho=0, zeta=-1, log=TRUE),
                 log(invprobit(xg, rho=0, zeta=-1)))

    # check oneminus
    expect_equal(invprobit(xg, rho=0, zeta=-1, oneminus = TRUE),
                 pnorm(xg, mean=0, sd=1, lower.tail = TRUE))
    expect_equal(invprobit(xg, rho=0, zeta=-1, oneminus = TRUE),
                 1 - invprobit(xg, rho=0, zeta=-1))
    # reducing the range because the manual calculation is too imprecise
    xg <- seq(-5, 5, length.out=101)
    expect_equal(invprobit(xg, rho=0, zeta=-1, oneminus = TRUE, log=TRUE),
                 log(1 - invprobit(xg, rho=0, zeta=-1)))
})

test_that("invprobit fails on bad input", {
    xg <- seq(-10, 10, length.out = 1001)
    expect_error(invprobit(xg, rho=1:3, zeta=0.1))
    expect_silent(pnorm(xg, mean=1:3, sd=0.1))

    expect_error(invprobit(xg, rho=1, zeta=c(1, -1)))
})
