

# test image reduction
test_that("image reduction works",{
  tol <- 0.5
  test_img <- matrix(
    c(1, 1, 2, 2,
      1, 1, 2, 2,
      3, 3, 4, 4,
      3, 3, 4, 4),
    nrow = 4L,
    byrow = TRUE)

test_same <- reduce_img_size(
  test_img,
  h = 3L,
  w = 4L)
testthat::expect(identical(test_img, test_same)) })


