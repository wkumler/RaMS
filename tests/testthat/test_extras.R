
test_that("pmppm works as expected", {
  expect_gt(min(pmppm(118.0865)), 118.08)
  expect_lt(max(pmppm(118.0865)), 118.1)

  expect_equal(min(pmppm(100, 2.5)), 100-0.00025)
  expect_equal(max(pmppm(100, 2.5)), 100+0.00025)

  expect_lt(pmppm(100)[1], pmppm(100)[2])

  expect_gt(max(pmppm(100, 10)), max(pmppm(100, 5)))

  expect_error(pmppm(100:110, 10))
})

test_that("pmppm returns a length-2 window around a single mass", {
  expect_equal(pmppm(100, 5), c(99.9995, 100.0005))
  expect_equal(pmppm(1000000, 5), c(999995, 1000005))
  expect_length(pmppm(118.0865), 2)
})

test_that("pmppm rejects vectorized input and points to pmppm_vec", {
  expect_error(pmppm(c(100, 200), 5), "pmppm_vec")
  expect_error(pmppm(100, c(5, 10)), "ppm")
  expect_error(pmppm("100", 5), "mass")
  expect_error(pmppm(100, "5"), "ppm")
  expect_warning(pmppm(100, -5), "positive")
})

test_that("pmppm_vec returns a named list of lower and upper bounds", {
  masses <- c(100, 118.0865, 1000000)
  win <- pmppm_vec(masses, 5)
  expect_type(win, "list")
  expect_named(win, c("lower", "upper"))
  expect_equal(win$lower, masses*(1-5/1e6))
  expect_equal(win$upper, masses*(1+5/1e6))
})

test_that("pmppm_vec matches pmppm for each mass", {
  masses <- c(90.0555, 104.1070, 118.0865, 892.535313)
  win <- pmppm_vec(masses, 2.5)
  for(i in seq_along(masses)){
    expect_identical(c(win$lower[i], win$upper[i]), pmppm(masses[i], 2.5))
  }
})

test_that("pmppm_vec returns a list for a single mass, with a message", {
  win <- pmppm_vec(100, 5)
  expect_type(win, "list")
  expect_identical(unlist(win, use.names=FALSE), pmppm(100, 5))
  expect_silent(pmppm_vec(c(100, 200), 5))
})

test_that("pmppm_vec handles empty and NA masses", {
  expect_identical(pmppm_vec(numeric(0), 5), list(lower=numeric(0), upper=numeric(0)))
  win <- pmppm_vec(c(100, NA), 5)
  expect_identical(is.na(win$lower), c(FALSE, TRUE))
  expect_identical(is.na(win$upper), c(FALSE, TRUE))
})

test_that("pmppm_vec uses the default ppm of 4", {
  expect_identical(pmppm_vec(c(100, 200)), pmppm_vec(c(100, 200), 4))
})

test_that("pmppm_vec validates its inputs", {
  expect_error(pmppm_vec(c("100", "200"), 5), "mass")
  expect_error(pmppm_vec(c(100, 200), "5"), "ppm")
  expect_error(pmppm_vec(c(100, 200), c(5, 10)), "ppm")
  expect_error(pmppm_vec(c(100, 200), NA_real_), "ppm")
  expect_error(pmppm_vec(c(100, 200), numeric(0)), "ppm")
  expect_warning(pmppm_vec(c(100, 200), -5), "positive")
})

test_that("pmppm_vec output works with data.table's %inrange%", {
  mz_vals <- c(90.0551, 104.1071, 104.2, 118.0866, 150)
  expect_identical(
    mz_vals %inrange% pmppm_vec(c(118.0865, 104.1070, 90.0555), 5),
    c(TRUE, TRUE, FALSE, TRUE, FALSE)
  )
  # Same answer as checking each window with pmppm and combining
  any_window <- Reduce(`|`, lapply(c(118.0865, 104.1070, 90.0555), function(m){
    mz_vals %between% pmppm(m, 5)
  }))
  expect_identical(mz_vals %inrange% pmppm_vec(c(118.0865, 104.1070, 90.0555), 5),
                   any_window)
})

test_that("pmppm_vec output works with data.table's %between% per element", {
  mz_vals <- c(90.0551, 104.1071, 104.2, 118.0866, 118.0866)
  targets <- c(90.0555, 104.1070, 104.1070, 118.0865, 90.0555)
  expect_identical(mz_vals %between% pmppm_vec(targets, 5),
                   c(TRUE, TRUE, FALSE, TRUE, FALSE))
  # Mismatched lengths error rather than recycling silently
  expect_error(mz_vals %between% pmppm_vec(c(100, 200), 5))
})

test_that("pmppm_vec bounds can be assigned as columns for a non-equi join", {
  compounds <- data.table::data.table(compound=c("betaine", "choline"),
                                      mass=c(118.0865, 104.1070))
  compounds[, c("lower", "upper"):=pmppm_vec(mass, 5)]
  expect_equal(compounds$lower, c(118.0865, 104.1070)*(1-5/1e6))
  expect_equal(compounds$upper, c(118.0865, 104.1070)*(1+5/1e6))

  pts <- data.table::data.table(mz=c(90.0551, 104.1071, 118.0866, 150))
  joined <- pts[compounds, .(compound, mz=x.mz), on=.(mz>=lower, mz<=upper),
                nomatch=NULL]
  expect_identical(joined$compound, c("betaine", "choline"))
  expect_identical(joined$mz, c(118.0866, 104.1071))
})

test_that("pmppm_vec works in the documented tidyverse patterns", {
  skip_if_not_installed("dplyr")
  ms1 <- data.frame(rt=1:6, mz=c(104.1071, 104.2, 118.0866, 118.0864, 136.0619, 150))
  cmpds <- data.frame(cmpd_name=c("betaine", "choline", "adenine"),
                      cmpd_mz=c(118.0865, 104.1070, 136.0618))

  filtered <- dplyr::filter(ms1, mz %inrange% pmppm_vec(cmpds$cmpd_mz, 5))
  expect_identical(filtered$rt, c(1L, 3L, 4L, 5L))

  windows <- dplyr::mutate(cmpds, as.data.frame(pmppm_vec(cmpd_mz, 5)))
  expect_named(windows, c("cmpd_name", "cmpd_mz", "lower", "upper"))
  joined <- dplyr::inner_join(windows, ms1,
                              dplyr::join_by(between(y$mz, x$lower, x$upper)))
  expect_identical(joined$cmpd_name, c("betaine", "betaine", "choline", "adenine"))
  expect_identical(joined$mz, c(118.0866, 118.0864, 104.1071, 136.0619))

  # One window per row comes from mapping pmppm over each mass (lapply here,
  # purrr's map() in the docs), and agrees with pmppm_vec
  per_row <- dplyr::mutate(cmpds, mz_range=lapply(cmpd_mz, pmppm, ppm=5))
  vec <- pmppm_vec(cmpds$cmpd_mz, 5)
  expect_identical(per_row$mz_range, Map(c, vec$lower, vec$upper))
  eic_rows <- vapply(per_row$mz_range, function(mzr)nrow(dplyr::filter(ms1, mz %between% mzr)), integer(1))
  expect_identical(eic_rows, c(2L, 1L, 1L))
})

test_that("pmppm_vec isn't a per-row list-column (documented pitfall)", {
  skip_if_not_installed("dplyr")
  three <- data.frame(mz=c(118.0865, 104.1070, 136.0618))
  expect_error(dplyr::mutate(three, mz_range=pmppm_vec(mz)), "size")
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
  test_df$startime <- rep(gl(2, 5, labels = c("Morn", "Eve")), each=100)

  expect_null(qplotMS1data(test_df, force_base = TRUE))
  expect_null(qplotMS1data(test_df, color_col="startime", force_base = TRUE))
  expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(118.0865)],
                               force_base = TRUE))
  expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(138.0555, 10)],
                               force_base = TRUE))

  if(requireNamespace("ggplot2", quietly = TRUE)){
    expect_s3_class(qplotMS1data(test_df), "ggplot")
    expect_no_error(print(qplotMS1data(test_df, color_col="startime") +
                            ggplot2::facet_wrap(~startime, ncol=2, scales="free")))
    expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(118.0865)]))
    expect_no_error(qplotMS1data(mzML_everything$MS1[mz%between%pmppm(138.0555, 10)]))
  }
})

# Two files, two m/z groups, with several points per scan
ms1_example <- data.frame(
  filename=rep(c("a.mzML", "b.mzML"), each=8),
  rt=rep(rep(1:2, each=4), 2),
  mz_group=rep(c("119.08", "119.08", "119.09", "119.09"), 4),
  int=c(10, 30, 5, 2,   20, 20, 1, 7,     # file a: rt 1, rt 2 (tie at 20)
        40, 15, 9, 12,  50, 60, 8, 3)     # file b: rt 1, rt 2
)
ms2_example <- data.frame(
  fragmz=c(58.0651, 59.0730, 60.0808, 118.0863, 58.0652, 118.0864),
  int=c(100, 15, 30, 60, 80, 90),
  rt=rep(c(1, 2), c(4, 2)),
  filename=rep(c("a.mzML", "b.mzML"), c(4, 2)),
  voltage=rep(c(20L, 35L), c(4, 2))
)

test_that("slice_max keeps the most intense point per scan and trace", {
  keep <- slice_max_int(ms1_example)
  expect_identical(ms1_example$int[keep], c(30, 20, 40, 60))
  # Rows 5 and 6 tie at 20; the first one wins
  expect_identical(which(keep)[2], 5L)
  # Discrete color columns split traces, continuous ones don't
  expect_identical(ms1_example$int[slice_max_int(ms1_example, "mz_group")],
                   c(30, 5, 20, 7, 40, 12, 60, 8))
  numeric_df <- transform(ms1_example, mz=seq(119.08, 119.09, length.out=16))
  expect_identical(slice_max_int(numeric_df, "mz"), keep)
})

test_that("slice_max matches dplyr::slice_max", {
  skip_if_not_installed("dplyr")
  ref <- dplyr::slice_max(ms1_example, int, by=c(rt, filename, mz_group), with_ties=FALSE)
  kept <- ms1_example[slice_max_int(ms1_example, "mz_group"),]
  expect_setequal(paste(kept$filename, kept$rt, kept$mz_group, kept$int),
                  paste(ref$filename, ref$rt, ref$mz_group, ref$int))
})

test_that("base plots draw with no, discrete, integer, and continuous colors", {
  pdf(NULL)
  on.exit(dev.off())
  depth_df <- transform(ms1_example, depth=ifelse(filename=="a.mzML", 5, 125))
  expect_null(qplotMS1data(ms1_example, force_base=TRUE))
  expect_null(qplotMS1data(ms1_example, color_col=NULL, force_base=TRUE))
  expect_null(qplotMS1data(ms1_example, color_col="mz_group", force_base=TRUE, slice_max=TRUE))
  expect_null(qplotMS1data(depth_df, color_col="depth", force_base=TRUE))
  expect_null(qplotMS1data(ms1_example[, c("rt", "int")], force_base=TRUE))
  expect_null(qplotMS2data(ms2_example, force_base=TRUE))
  expect_null(qplotMS2data(ms2_example, color_col="voltage", force_base=TRUE, relative=TRUE))
  expect_error(qplotMS1data(ms1_example, color_col="nope"), "nope")
})

test_that("ggplot draws one line per file and discrete color level", {
  skip_if_not_installed("ggplot2")
  n_groups <- function(gg)length(unique(ggplot2::ggplot_build(gg)$data[[1]]$group))
  expect_identical(n_groups(qplotMS1data(ms1_example)), 2L)
  expect_identical(n_groups(qplotMS1data(ms1_example, color_col="mz_group")), 4L)
  # Integers are discrete, doubles continuous
  int_df <- transform(ms1_example, mz_group=as.integer(factor(mz_group)))
  expect_identical(n_groups(qplotMS1data(int_df, color_col="mz_group")), 4L)
  depth_df <- transform(ms1_example, depth=ifelse(filename=="a.mzML", 5, 125))
  expect_identical(n_groups(qplotMS1data(depth_df, color_col="depth")), 2L)
})

test_that("qplotMS2data draws stems to zero and relative=TRUE scales each scan", {
  skip_if_not_installed("ggplot2")
  built <- ggplot2::ggplot_build(qplotMS2data(ms2_example))
  expect_identical(built$data[[1]]$x, built$data[[1]]$xend)
  expect_true(all(built$data[[1]]$yend==0))
  rel <- qplotMS2data(ms2_example, relative=TRUE)
  expect_equal(ggplot2::ggplot_build(rel)$data[[2]]$y, c(100, 15, 30, 60, 80/90*100, 100))
  expect_identical(rel$labels$y, "Relative intensity (%)")
  expect_identical(ms2_example$int[1:2], c(100, 15))
})

test_that("color_col only defaults to filename when there are multiple files", {
  skip_if_not_installed("ggplot2")
  has_color <- function(gg)"colour"%in%names(gg$mapping)
  expect_true(has_color(qplotMS1data(ms1_example)))
  expect_false(has_color(qplotMS1data(ms1_example[ms1_example$filename=="a.mzML",])))
  expect_false(has_color(qplotMS2data(ms2_example[ms2_example$filename=="a.mzML",])))
  expect_true(has_color(qplotMS2data(ms2_example[ms2_example$filename=="a.mzML",], color_col="filename")))
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
