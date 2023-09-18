## Stores all functions used to complete analysis that are not in the SEraster package.

calculateDensity <- function(matrix.array) {
  sum(matrix.array != 0)/(dim(matrix.array)[1] * dim(matrix.array)[2])
}