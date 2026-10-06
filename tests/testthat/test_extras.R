
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

  baseplot_output <- qplotMS1data(test_df, force_base = TRUE)
  expect_null(baseplot_output)

  if(requireNamespace("ggplot2", quietly = TRUE)){
    ggplot_output <- qplotMS1data(test_df)
    expect_s3_class(ggplot_output, "ggplot")


    test_df$startime <- rep(gl(2, 5, labels = c("Morn", "Eve")), each=100)
    expect_no_error(qplotMS1data(test_df, color_col="startime", facet_col="startime"))
    expect_no_error({
      qplotMS1data(test_df, color_col="startime", facet_col="startime",
                   facet_args=list(ncol=2, scales="free"))
    })

    expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(118.0865)]))
    expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(138.0555, 10)]))
  }
  expect_warning(qplotMS1data(test_df, color_col="startime", force_base = TRUE))
  expect_warning(qplotMS1data(test_df, facet_col="startime", force_base = TRUE))

  expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(118.0865)],
                               force_base = TRUE))
  expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(138.0555, 10)],
                               force_base=TRUE))
})

example_mz_vals <- c(118.0, 118.1, 138.0, 152.0, 118.2, 138.1, 118.1)

test_that("mz_group works as expected", {
  expect_equal(mz_group(example_mz_vals, ppm = 1), c(1, 2, 3, 4, 5, 6, 2))
  expect_equal(mz_group(example_mz_vals, ppm = 1000), c(1, 1, 2, 3, 4, 2, 1))
  expect_equal(mz_group(example_mz_vals, ppm = 2e5), c(1, 1, 1, 2, 1, 1, 1))

  expect_equal(mz_group(example_mz_vals, ppm = 1000, min_group_size = 2),
               c(1, 1, 2, NA, NA, 2, 1))
  expect_equal(mz_group(example_mz_vals, ppm = 1000, max_groups = 2),
               c(1, 1, 2, NA, NA, 2, 1))
  expect_equal(max(mz_group(runif(1000), ppm = 1e5, max_groups=2), na.rm = TRUE), 2)

  mz_vals <- mzML_everything$MS1[mz%between%pmppm(119.0865, 100)][order(int, decreasing = TRUE)]$mz
  expect_equal(max(mz_group(mz_vals, ppm=10)), 6)
  expect_equal(max(mz_group(mz_vals, ppm=50)), 2)
  expect_equal(max(mz_group(mz_vals, ppm=5, max_groups = 2), na.rm = TRUE), 2)
  expect_equal(max(mz_group(mz_vals, ppm=5, max_groups = 3), na.rm = TRUE), 3)
  expect_equal(max(mz_group(mz_vals, ppm=5, min_group_size = 1), na.rm = TRUE), 11)
  expect_true(any(is.na(mz_group(mz_vals, ppm=5, max_groups = 2))))
  # expect_true(any(is.na(mz_group(mz_vals, ppm=5, min_group_size = 1))))
})

# Simulated MS1-like data: a few hundred true masses with ~1.5 ppm scatter
# plus uniform noise, ordered by random "intensity"
sim_mz <- function(n_masses=300, pts_per=50, noise_n=3000, seed=1){
  set.seed(seed)
  m <- runif(n_masses, 70, 1000)
  mz <- c(rep(m, each=pts_per)*(1+rnorm(n_masses*pts_per, 0, 1.5e-6)),
          runif(noise_n, 70, 1000))
  mz[order(rlnorm(length(mz), 10, 2), decreasing=TRUE)]
}

test_that("mz_group returns integer groups for the documented examples", {
  expect_identical(mz_group(example_mz_vals, ppm = 1), c(1L, 2L, 3L, 4L, 5L, 6L, 2L))
  expect_identical(mz_group(example_mz_vals, ppm = 1000), c(1L, 1L, 2L, 3L, 4L, 2L, 1L))
  # 20% of 118 reaches 141.6, so 152 is still on its own
  expect_identical(mz_group(example_mz_vals, ppm = 200000), c(1L, 1L, 1L, 2L, 1L, 1L, 1L))
})

test_that("mz_group uses the first unassigned value as each window center", {
  # Centered on 118.0, 118.2 falls outside 1000 ppm (+/-0.118)...
  expect_identical(mz_group(c(118.0, 118.1, 118.2), 1000), c(1L, 1L, 2L))
  # ...but centered on 118.1 everything is within +/-0.1181
  expect_identical(mz_group(c(118.1, 118.0, 118.2), 1000), c(1L, 1L, 1L))
})

test_that("mz_group windows are inclusive at the ppm boundary", {
  edge <- 100*(1+5/1e6)
  expect_identical(mz_group(c(100, edge), ppm = 5), c(1L, 1L))
  expect_identical(mz_group(c(100, edge*(1+1e-9)), ppm = 5), c(1L, 2L))
})

test_that("mz_group treats min_group_size as an inclusive minimum", {
  # Window sizes for example_mz_vals at 1000 ppm are 3, 2, 1, 1
  expect_identical(mz_group(example_mz_vals, 1000, min_group_size = 0), mz_group(example_mz_vals, 1000))
  expect_identical(mz_group(example_mz_vals, 1000, min_group_size = 1), mz_group(example_mz_vals, 1000))
  expect_identical(mz_group(example_mz_vals, 1000, min_group_size = 2),
                   c(1L, 1L, 2L, NA, NA, 2L, 1L))
  expect_identical(mz_group(example_mz_vals, 1000, min_group_size = 3),
                   c(1L, 1L, NA, NA, NA, NA, 1L))
  expect_identical(mz_group(example_mz_vals, 1000, min_group_size = 4), rep(NA_integer_, 7))
})

test_that("mz_group numbers kept groups consecutively", {
  x <- sim_mz()
  for(mgs in c(1, 3, 20)){
    g <- mz_group(x, 5, min_group_size = mgs)
    expect_identical(sort(unique(g[!is.na(g)])), seq_len(max(g, na.rm = TRUE)))
    expect_true(all(tabulate(g) >= mgs))
  }
})

test_that("mz_group max_groups limits only kept groups", {
  expect_identical(mz_group(example_mz_vals, 1000, max_groups = 1),
                   c(1L, 1L, NA, NA, NA, NA, 1L))
  expect_identical(mz_group(example_mz_vals, 1000, max_groups = 0), rep(NA_integer_, 7))
  # A singleton found first is dropped and doesn't use up the quota
  expect_identical(mz_group(c(152, 118, 118.1, 138, 138.1), 1000,
                            min_group_size = 2, max_groups = 2),
                   c(NA, 1L, 1L, 2L, 2L))
})

test_that("mz_group returns all groups when max_groups exceeds the number available", {
  expect_identical(mz_group(example_mz_vals, 1000, max_groups = 100), mz_group(example_mz_vals, 1000))
  expect_identical(mz_group(example_mz_vals, 1000, min_group_size = 3, max_groups = 5),
                   c(1L, 1L, NA, NA, NA, NA, 1L))
})

test_that("mz_group max_groups equals the unlimited result truncated", {
  x <- sim_mz()
  full <- mz_group(x, 5, min_group_size = 3)
  for(k in c(1, 7, 50)){
    expected <- full
    expected[expected > k] <- NA
    expect_identical(mz_group(x, 5, min_group_size = 3, max_groups = k), expected)
  }
})

test_that("mz_group handles NA and empty input", {
  expect_identical(mz_group(c(118, NA, 118.0001), 5), c(1L, NA, 1L))
  expect_identical(mz_group(c(NA_real_, NA_real_), 5), c(NA_integer_, NA_integer_))
  expect_identical(mz_group(numeric(0), 5), integer(0))
})

test_that("mz_group validates its inputs", {
  expect_error(mz_group(as.character(example_mz_vals), 5), "mz_vals")
  expect_error(mz_group(example_mz_vals, "5"), "ppm")
  expect_error(mz_group(example_mz_vals, c(5, 10)), "ppm")
  expect_error(mz_group(example_mz_vals, 0), "ppm")
  expect_error(mz_group(example_mz_vals, 5, min_group_size = "2"), "min_group_size")
  expect_error(mz_group(example_mz_vals, 5, max_groups = -1), "max_groups")
  expect_error(mz_group(example_mz_vals, 5, max_groups = c(1, 2)), "max_groups")
  expect_error(mz_group(example_mz_vals, 5, mz_factor = NA), "mz_factor")
  expect_error(mz_group(example_mz_vals, 5, mz_factor = "yes"), "mz_factor")
})

test_that("mz_group mz_factor labels groups by mean m/z, levels sorted by mass", {
  f <- mz_group(example_mz_vals, 1000, mz_factor = TRUE)
  expect_s3_class(f, "factor")
  expect_identical(levels(f), c("118.06667", "118.20000", "138.05000", "152.00000"))
  expect_identical(as.character(f),
                   c("118.06667", "118.06667", "138.05000", "152.00000",
                     "118.20000", "138.05000", "118.06667"))
})

test_that("mz_group mz_factor matches the integer grouping", {
  x <- sim_mz()
  g <- mz_group(x, 5, min_group_size = 3)
  f <- mz_group(x, 5, min_group_size = 3, mz_factor = TRUE)
  expect_identical(is.na(f), is.na(g))
  expect_identical(nlevels(f), max(g, na.rm = TRUE))
  # One-to-one mapping between integer groups and factor levels
  expect_true(all(rowSums(table(g, f) > 0) == 1))
  expect_identical(as.character(f[!is.na(g)]),
                   sprintf("%.5f", ave(x, g, FUN = mean))[!is.na(g)])
  expect_false(is.unsorted(as.numeric(levels(f))))
})

test_that("mz_group mz_factor respects min_group_size and max_groups", {
  expect_identical(as.character(mz_group(example_mz_vals, 1000, min_group_size = 3, mz_factor = TRUE)),
                   c("118.06667", "118.06667", NA, NA, NA, NA, "118.06667"))
  expect_identical(levels(mz_group(example_mz_vals, 1000, max_groups = 2, mz_factor = TRUE)),
                   c("118.06667", "138.05000"))
})

test_that("mz_group mz_factor handles NA, empty, and no kept groups", {
  expect_identical(mz_group(numeric(0), 5, mz_factor = TRUE), factor(character(0)))
  all_na <- mz_group(c(NA_real_, NA_real_), 5, mz_factor = TRUE)
  expect_s3_class(all_na, "factor")
  expect_true(all(is.na(all_na)))
  none <- mz_group(example_mz_vals, 1, min_group_size = 10, mz_factor = TRUE)
  expect_s3_class(none, "factor")
  expect_true(all(is.na(none)))
  expect_identical(nlevels(none), 0L)
})

test_that("mz_group mz_factor merges groups whose means round to the same label", {
  f <- mz_group(c(100, 100.000001), ppm = 0.001, mz_factor = TRUE)
  expect_identical(mz_group(c(100, 100.000001), ppm = 0.001), c(1L, 2L))
  expect_identical(levels(f), "100.00000")
  expect_identical(as.integer(f), c(1L, 1L))
})
