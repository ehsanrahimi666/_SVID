test_that("algebraically identical expressions share a canonical form", {
  e <- c("(N / B) / R", "(N / R) / B", "N / (B * R)")
  cf <- svid_canonicalize(e)
  expect_length(unique(cf), 1L)
})

test_that("distinct rational functions do not collide", {
  cf <- svid_canonicalize(c("N / R", "R / N", "(N - R) / (N + R)", "N - R"))
  expect_length(unique(cf), 4L)
})

test_that("constant multiples are algebraically distinct", {
  expect_false(identical(svid_canonicalize("2 * (N / G)"),
                         svid_canonicalize("N / G")))
})

test_that("nested and flattened forms agree", {
  expect_identical(svid_canonicalize("((N - R)) / ((N + R))"),
                   svid_canonicalize("(N - R) / (N + R)"))
})

test_that("unsupported operators are reported as NA rather than silently wrong", {
  expect_true(is.na(svid_canonicalize("log(N / R)")))
})
