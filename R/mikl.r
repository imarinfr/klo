# Gaussian estimate of the differential entropy of the multivariate
# random variable x. If x is a multivariate Gaussian, then the estimate is
# asymptotically unbiased, otherwise it is just an approximation.
#
# x is a n-by-d numeric matrix, in which the n rows correspond to
# observations and the d columns to variables (or coordinates) of the
# multivariate distributions
entg <- function(x) {
  # convert x to matrix if not already
  x <- as.matrix(x)
  # get dimensions
  n <- nrow(x)
  d <- ncol(x)
  if(qr(x)$rank != d) stop("covariance matrix of x does not have full rank")
  return(1 / 2 * log2((2 * pi * exp(1))^d * det(cov(x))))
}

# Gaussian estimate of mutual information between the multivariate
# random variables x and y. If both x and y are multivariate Gaussians,
# then the estimate is asymptotically unbiased, otherwise it is just an
# approximation
#
# x, y are each n-by-d numeric matrices, in which the n rows correspond to
# observations and the d columns to variables (or coordinates) of the
# multivariate distributions
mig <- function(x, y)
  # computation of mutual information is based on its decomposition in
  # marginal and joint differential entropies
  return(entg(x) + entg(y) - entg(cbind(x, y)))

# The Kozachenko-Leonenko (KL) estimate of the differential entropy
# of the multivariate random variable x.
#
# If type is 'kl', the raw KL estimator is computed.
# If type is 'klo', the offset KL estimator is computed.
#
# k is the nearest neighbor for which to search, 1 the nearest neighbor, 2,
# the second nearest, and so on
#
# x is a n-by-d numeric matrix, in which the n rows correspond to
# observations and the d columns to variables (or coordinates) of the
# multivariate distributions
entkl <- function(x, type = "klo", k = 1) {
  require(FNN) # need a fast nearest-neighbor algorithm
  require(expm)
  # check input
  if(type != "kl" && type != "klo")
    stop("incorrect type of estimator")
  if(k <= 0)
    stop("k-nearest neighbour 'k' must be larger than 0")
  # convert x to matrix if not already
  x <- as.matrix(x)
  # get dimensions
  n <- nrow(x)
  d <- ncol(x)
  if(qr(x)$rank != d) stop("covariance matrix of x does not have full rank")
  # if type is 'kl' then calculate differential entropy with x directly.
  # If type is 'klo', then calculate differential entropy as if x were
  # Gaussian and prepare the data to calculate the offset differential
  # entropy using nearest-neighbor algorithms
  entgp <- 0
  if(type == "klo") {
    entgp <- entg(x)
    x     <- x %*% solve(sqrtm(cov(x))) / sqrt(2 * pi * exp(1))
  }
  # the nearest-neighbour distances
  ldist <- log(get.knnx(x, x, k = k + 1, algorithm ="kd_tree")$nn.dist[,k + 1])
  ldist[ldist == -Inf] <- 0
  # differential entropy calculation in nats
  ent <- d * mean(ldist) + log(n - 1) - psigamma(k) + log(2 * pi^(d / 2) / (d * gamma(d / 2)))
  # return differential entropy in bits
  return(log2(exp(ent)) + entgp)
}

# The Kozachenko-Leonenko (KL) estimate of mutual information between
# the multivariate random variables x and y.
#
# If type is 'kl', the raw KL estimator is computed.
# If type is 'klo', the offset KL estimator is computed.
#
# k is the nearest neighbor for which to search, 1 the nearest neighbor, 2,
# the second nearest, and so on
#
# x, y are each n-by-d numeric matrices, in which the n rows correspond to
# observations and the d columns to variables (or coordinates) of the
# multivariate distributions
mikl <- function(x, y, type = "klo", k = 1)
  # computation of mutual information is based on its decomposition in
  # marginal and joint differential entropies
  return(entkl(x, type, k) + entkl(y, type, k) - entkl(cbind(x, y), type, k))