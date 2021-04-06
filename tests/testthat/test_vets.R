context("Tests for ves() function")

Y <- cbind(Mcomp::M3$N2570$x,Mcomp::M3$N2571$x)

# Basic VES check
testModel <- suppressWarnings(vets(Y,"MMdM", silent=TRUE))
test_that("Test VETS(MMdM)", {
    expect_match(testModel$model, "MMdM")
})

# Test VES with everything individual
test_that("Test VES with all individual", {
    skip_on_cran()
    testModel <- vets(Y,"MMdM", parameters="none", initials="none", components="none", silent=TRUE)
    expect_equal(length(coefficients(testModel)), 34)
})

# Test VES with common seasonal and persistence
test_that("Test VETS with common seasonal and persistence", {
    skip_on_cran()
    testModel <- vets(Y,"MMdM", components=c("seasonal"), parameters=c("l","t","s"), silent=TRUE)
    expect_equal(length(coefficients(testModel)), 20)
})

# Test VETS with common initials
test_that("Test VETS with common initials", {
    skip_on_cran()
    testModel <- vets(Y,"AAN", initials=c("l","t"), parameters="none", silent=TRUE)
    expect_equal(length(coefficients(testModel)), 6)
})

# Test VETS with a trace cost function
test_that("Test VETS with a trace cost function", {
    skip_on_cran()
    testModel <- vets(Y,"AAN", loss="t", silent=TRUE)
    expect_match(testModel$loss, "trace")
})

# Model selection in VETS
test_that("Model selection in VETS", {
    skip_on_cran()
    testModel <- vets(Y, "PPP", silent=TRUE)
    expect_match(testModel$loss, "likelihood")
})