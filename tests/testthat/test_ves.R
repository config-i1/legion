context("Tests for ves() function")

Y <- ts(cbind(1000+0.5*c(1:100)+rnorm(100,0,10),
              cbind(1000+1.5*c(1:100)+rnorm(100,0,10))),
        frequency=12)

# Basic VES check
testModel <- suppressWarnings(ves(Y,"MMdM", silent=TRUE))
test_that("Test VES(MMdM)", {
    expect_match(testModel$model, "MMdM")
})

# Reuse previous VES
test_that("Reuse VES", {
    expect_equal(ves(Y, model=testModel, silent=TRUE)$persistence, testModel$persistence)
})

# Test VES with individual seasonality and persistence
test_that("Test VES with individual seasonality and persistence", {
    skip_on_cran()
    testModel <- ves(Y,"MMdM", initialSeason="i", persistence="i", silent=TRUE)
    expect_equal(length(coefficients(testModel)), 35)
})

# Test VES with grouped initials and dependent persistence
test_that("Test VES with grouped initials and dependent persistence", {
    skip_on_cran()
    testModel <- ves(Y,"AAN", initial="c", persistence="d", silent=TRUE)
    expect_equal(length(coefficients(testModel)), 10)
})

# Test VES with a trace cost function
test_that("Test VES with a trace cost function", {
    skip_on_cran()
    testModel <- ves(Y,"AAN", loss="t", silent=TRUE)
    expect_match(testModel$loss, "trace")
})

# Test VES with a dependent transition and independent interval
test_that("Test VES with a dependent transition and independent interval", {
    skip_on_cran()
    testModel <- ves(Y,"AAN", transition="d", silent=TRUE)
    expect_false(isTRUE(all.equal(testModel$transition[1,4], 0)))
    testForecast <- forecast(testModel, h=10, interval="prediction")
    expect_equal(dim(testForecast$PI),c(10,4))
})

# Model selection in VES
test_that("Model selection in VES", {
    skip_on_cran()
    testModel <- ves(Y,"PPP", silent=TRUE)
    expect_match(testModel$loss, "likelihood")
})

# Simulate the data from VES
testModel <- ves(Y,"AAN", silent=TRUE)
test_that("VES based on pre-estimated model", {
    expect_match(substr(simulate(testModel,nsim=10,obs=100)$model,1,3), "VES")
})
