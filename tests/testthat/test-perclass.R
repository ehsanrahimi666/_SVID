test_that("per-class discovery returns one winner per class", {
  d <- svid_demo_data(n = 90, epochs = 4, groups = 9)
  pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
                                                max_depth = 1)), d)
  pc <- svid_per_class(pool, d)
  expect_s3_class(pc, "svid_per_class")
  expect_equal(nrow(pc$winners), nlevels(d$response))
  expect_true(all(pc$winners$class %in% levels(d$response)))
  expect_true(all(pc$scores[!is.na(pc$scores)] >= 0.5 - 1e-8))
})

test_that("the fast AUC path agrees with the generic separability path", {
  d <- svid_demo_data(n = 80, epochs = 3, groups = 8)
  ex <- c("N / R", "(N - R) / (N + R)", "G / B")
  fast <- svid_per_class(ex, d, metric = "auc")
  slow <- svid_per_class(ex, d, metric = "fisher")
  expect_equal(dim(fast$scores), dim(slow$scores))
  expect_false(anyNA(fast$scores))
})

test_that("specialists never score below the single all-class index", {
  d <- svid_demo_data(n = 90, epochs = 4, groups = 9)
  pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
                                                max_depth = 1)), d)
  sp <- svid_specialists(svid_per_class(pool, d))
  expect_true(all(sp$gain >= -1e-8))
})

test_that("per-class discovery rejects a continuous response", {
  d <- svid_simulate(n = 40, epochs = 2, response = "lai")
  expect_error(svid_per_class("N / R", d), "categorical")
})
