
test_that("pmppm works as expected", {
  expect_gt(min(pmppm(118.0865)), 118.08)
  expect_lt(max(pmppm(118.0865)), 118.1)

  expect_equal(min(pmppm(100, 2.5)), 100-0.00025)
  expect_equal(max(pmppm(100, 2.5)), 100+0.00025)

  expect_lt(pmppm(100)[1], pmppm(100)[2])

  expect_gt(max(pmppm(100, 10)), max(pmppm(100, 5)))

  expect_warning(pmppm(100:110, 10))
})

test_that("trapz works as expected", {
  expect_equal(trapz(1:10, 1:10), 49.5)
  expect_equal(trapz(1:10, 10:1), 49.5)
  expect_equal(trapz(0, 1), 0)
  expect_equal(trapz(1, 1), 0)
  expect_equal(trapz(1:10, 1:10, baseline = "square"), 40.5)
  expect_equal(trapz(1:10, 1:10, baseline = "trapezoid"), 0)
  expect_equal(trapz(c(0, 100), c(1, 1)), 100)
  expect_equal(trapz(c(0, 100), c(1, 1), "square"), 0)
  expect_equal(trapz(c(0, 100), c(1, 1), "trapezoid"), 0)

  x_vals <- 1:10
  y_vals <- runif(10)
  expect_equal(trapz(x_vals, y_vals), trapz(x_vals, rev(y_vals)))
  expect_false(trapz(x_vals, y_vals)==trapz(x_vals, sample(y_vals)))

  expect_error(trapz(10))
  expect_error(trapz(1, 1:10))
  expect_warning(trapz(c(1, 1:10), c(1:11)))
  expect_warning(trapz(c(11, 1:10), c(1:11)))
  expect_error(trapz("banana", "apple"))
  expect_error(trapz(1:10, 1:10, "apple"))
  expect_warning(trapz(1:10000, runif(10000)))
  expect_warning(trapz(1:11, -5:5))
})
