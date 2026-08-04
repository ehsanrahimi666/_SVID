test_that("monotone transformations fall into one class", {
  d <- svid_demo_data(n = 60, epochs = 3)
  eq <- svid_equivalence(c("N / R", "R / N", "(N - R) / (N + R)", "N - R"), d)
  expect_equal(eq$n_classes, 2L)
  expect_true(all(c("N / R", "R / N", "(N - R) / (N + R)") %in%
                    eq$members[["N / R"]]))
})

test_that("invariance = none disables the reduction", {
  d <- svid_demo_data(n = 40, epochs = 2)
  eq <- svid_equivalence(c("N / R", "R / N"), d, invariance = "none")
  expect_equal(eq$n_classes, 2L)
})

test_that("pruning records a monotone reduction history", {
  d <- svid_demo_data(n = 60, epochs = 3)
  p <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"))), d)
  h <- svid_redundancy(p)
  expect_true(all(diff(h$n) <= 0))
  expect_equal(length(p), h$n[nrow(h)])
})
