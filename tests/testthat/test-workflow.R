test_that("generation respects the grammar", {
  g <- svid_grammar(c("B", "G", "R", "N"), ops = c("/"), max_depth = 1)
  p <- svid_generate(g)
  expect_true(all(grepl("/", as.character(p))))
  expect_false(any(grepl("\\+", as.character(p))))
})

test_that("oversized grammars are refused with a helpful message", {
  g <- svid_grammar(c("B", "G", "R", "N", "S1", "S2"), max_depth = 4)
  expect_error(svid_generate(g, max_size = 1000), "grammar would generate")
})

test_that("evaluation returns one matrix per expression", {
  d <- svid_demo_data(n = 40, epochs = 4)
  v <- svid_eval(c("N / R", "N - R"), d)
  expect_length(v, 2L)
  expect_equal(dim(v[[1]]), c(40L, 4L))
})

test_that("rank cleaning keeps values finite and order-preserving", {
  x <- matrix(c(1, Inf, -Inf, NA, 2), 5, 1)
  y <- svid:::.clean_feature(x, "rank")
  expect_true(all(is.finite(y)))
  expect_gt(y[2], y[5])
})

test_that("grouped resampling never splits a group", {
  d <- svid_demo_data(n = 60, groups = 6)
  rs <- svid_resample(d, "group", folds = 3)
  for (s in rs) {
    expect_length(intersect(d$meta$group[s$train], d$meta$group[s$test]), 0L)
  }
})

test_that("random splitting of grouped data warns", {
  d <- svid_demo_data(n = 60, groups = 6)
  expect_warning(svid_resample(d, "random"), "group")
})

test_that("evaluation produces metrics and predictions", {
  d <- svid_demo_data(n = 60, epochs = 3, groups = 6)
  fit <- svid_evaluate(c("(N - R) / (N + R)", "N / R"), d, n_seeds = 1,
                       baselines = "none")
  expect_s3_class(fit, "svid_fit")
  expect_true(all(c("oa", "kappa", "macro_f1") %in% names(fit$metrics)))
  expect_true(all(fit$metrics$oa >= 0 & fit$metrics$oa <= 100))
})

test_that("regression targets are handled", {
  d <- svid_simulate(n = 60, epochs = 3, response = "lai", groups = 6)
  expect_equal(d$type, "regression")
  fit <- svid_regress(c("(N - R) / (N + R)"), d, n_seeds = 1, baselines = "none")
  expect_true("r2" %in% names(fit$metrics))
})

test_that("confidence intervals bracket the estimate", {
  d <- svid_demo_data(n = 60, epochs = 3, groups = 6)
  fit <- svid_evaluate(c("(N - R) / (N + R)", "N / R"), d, n_seeds = 1,
                       baselines = "none")
  ci <- svid_ci(fit)
  expect_true(all(ci$lower <= ci$estimate & ci$estimate <= ci$upper))
})

test_that("cross-sensor portability detects missing bands", {
  x <- svid_cross_sensor("(RE2 - R) / (RE2 + R)", "sentinel2", "landsat8")
  expect_false(x$portable)
  expect_equal(x$missing_bands, "RE2")
})

test_that("the catalogue is filtered by available bands", {
  a <- svid_catalog(c("B", "G", "R", "N"))
  expect_true(all(vapply(strsplit(a$band_list, ","),
                         function(b) all(b %in% c("B", "G", "R", "N")),
                         logical(1))))
  expect_true("NDVI" %in% a$short_name)
})

test_that("temporal summaries have the requested columns", {
  d <- svid_demo_data(n = 30, epochs = 5)
  s <- svid_temporal(svid_eval("N / R", d)[[1]], stats = c("mean", "amplitude"))
  expect_equal(colnames(s), c("mean", "amplitude"))
})
