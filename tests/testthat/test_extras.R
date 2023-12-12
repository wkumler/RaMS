
test_that("pmppm works as expected", {
  expect_gt(min(pmppm(118.0865)), 118.08)
  expect_lt(max(pmppm(118.0865)), 118.1)

  expect_equal(min(pmppm(100, 2.5)), 100-0.00025)
  expect_equal(max(pmppm(100, 2.5)), 100+0.00025)

  expect_lt(pmppm(100)[1], pmppm(100)[2])

  expect_gt(max(pmppm(100, 10)), max(pmppm(100, 5)))

  expect_error(pmppm(100:110, 10))
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
  # expect_false(trapz(x_vals, y_vals)==trapz(x_vals, sample(y_vals)))

  expect_error(trapz(10))
  expect_error(trapz(1, 1:10))
  expect_warning(trapz(c(1, 1:10), c(1:11)))
  expect_warning(trapz(c(11, 1:10), c(1:11)))
  expect_error(trapz("banana", "apple"))
  expect_error(trapz(1:10, 1:10, "apple"))
  expect_warning(trapz(1:10000, runif(10000)))
  expect_warning(trapz(1:11, -5:5))
})

test_that("qplotMS1data works as expected", {
  test_df <- expand.grid(rt=rep(1:100, length.out=1000))
  test_df$int <- rep(dnorm(seq(-10, 10, length.out=100)), 10)*10+runif(1000)
  test_df$filename <- rep(LETTERS[1:10], each=100)

  ggplot_output <- qplotMS1data(test_df)
  expect_s3_class(ggplot_output, "ggplot")

  baseplot_output <- qplotMS1data(test_df, force_base = TRUE)
  expect_null(baseplot_output)

  test_df$startime <- rep(gl(2, 5, labels = c("Morn", "Eve")), each=100)
  expect_no_error(qplotMS1data(test_df, color_col="startime", facet_col="startime"))
  expect_no_error({
    qplotMS1data(test_df, color_col="startime", facet_col="startime",
                 facet_args=list(ncol=2, scales="free"))
  })
  expect_warning(qplotMS1data(test_df, color_col="startime", force_base = TRUE))
  expect_warning(qplotMS1data(test_df, facet_col="startime", force_base = TRUE))

  expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(118.0865)]))
  expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(138.0555, 10)]))
})

test_that("mz_group works as expected", {
  example_mz_vals <- c(118.0, 118.1, 138.0, 152.0, 118.2, 138.1, 118.1)
  expect_equal(mz_group(example_mz_vals, ppm = 1), c(1, 2, 3, 4, 5, 6, 2))
  expect_equal(mz_group(example_mz_vals, ppm = 1000), c(1, 1, 2, 3, 4, 2, 1))
  expect_equal(mz_group(example_mz_vals, ppm = 2e5), c(1, 1, 1, 2, 1, 1, 1))

  expect_equal(mz_group(example_mz_vals, ppm = 1000, min_group_size = 2),
               c(1, 1, NA, NA, NA, NA, 1))
  expect_equal(mz_group(example_mz_vals, ppm = 1000, max_groups = 2),
               c(1, 1, 2, NA, NA, 2, 1))
  expect_equal(max(mz_group(runif(1000), ppm = 1e5, max_groups=2), na.rm = TRUE), 2)

  mz_vals <- mzML_everything$MS1[mz%between%pmppm(119.0865, 100)][order(int, decreasing = TRUE)]$mz
  expect_equal(max(mz_group(mz_vals, ppm=10)), 5)
  expect_equal(max(mz_group(mz_vals, ppm=50)), 2)
  expect_equal(max(mz_group(mz_vals, ppm=5, max_groups = 2), na.rm = TRUE), 2)
  expect_equal(max(mz_group(mz_vals, ppm=5, max_groups = 3), na.rm = TRUE), 3)
  expect_equal(max(mz_group(mz_vals, ppm=5, min_group_size = 1), na.rm = TRUE), 4)
  expect_true(any(is.na(mz_group(mz_vals, ppm=5, max_groups = 2))))
  expect_true(any(is.na(mz_group(mz_vals, ppm=5, min_group_size = 1))))
})
