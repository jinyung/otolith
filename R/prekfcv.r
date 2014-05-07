# Generate folds
# @description Give the folds of k-fold cross-validation
# @author Matthias C. M. Troffaes
# @references https://gist.github.com/mcmtroffaes/709908

kfcv.sizes = function(n, k=10) {
  sizes = c()
  for (i in 1:k) {
    first = 1 + (((i - 1) * n) %/% k)
    last = ((i * n) %/% k)
    sizes = append(sizes, last - first + 1)
  }
  sizes
}

kfcv.testing = function(n, k=10) {
  indices = list()
  sizes = kfcv.sizes(n, k=k)
  values = 1:n
  for (i in 1:k) {
    s = sample(values, sizes[i])
    indices[[i]] = s
    values = setdiff(values, s)
  }
  indices
}