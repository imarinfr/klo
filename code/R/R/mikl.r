#' @importFrom stats cov
#' @importFrom expm sqrtm
#' @importFrom FNN get.knnx

#' @rdname mikl
#' @title Estimates of the differential entropy and mutual information
#' @description Functions that compute the Gaussian and the regular
#' and offset versions of the Kozachenko-Leonenko estimators of the
#' differential entropy and mutual information for multivariate
#' continuous variables
#' @section Differential entropy estimates:
#' \itemize{
#'   \item\code{entG}  Gaussian estimate of the differential entropy of the
#'                     multivariate random variables \code{x} and \code{y}.
#'                     If both \code{x} and \code{y} are multivariate Gaussians,
#'                     then the estimate is asymptotically unbiased, otherwise
#'                     it is just an approximation
#'   \item\code{entkl} The Kozachenko-Leonenko (KL) estimate of the
#'                     differential entropy of the multivariate random variable
#'                     \code{x}
#' }
#' @section Mutual information entropy estimates:
#' \itemize{
#'   \item\code{miG}   Gaussian estimate of mutual information between the
#'                     multivariate random variables \code{x} and \code{y}.
#'                     If both \code{x} and \code{y} are multivariate Gaussians,
#'                     then the estimate is asymptotically unbiased, otherwise
#'                     it is just an approximation
#'   \item\code{mikl}  The Kozachenko-Leonenko (KL) estimate of the
#'                     mutual information of the multivariate random variable
#'                     \code{x}
#' }
#' @param x,y n-by-d numeric matrices, in which the n rows correspond to observations
#'            and the d columns to variables (or coordinates) of the multivariate
#'            distributions
#' @param type is the type of estimator, \code{"kl"} for the Kozachenko-Leonenko and
#'             \code{"klo"} (default) for its offset version.
#' @param k is the rank of the nearest neighbor for which to search, \code{1},
#'          the nearest neighbor (default), \code{2}, the second nearest,
#'          and so on.
#' @examples
#' # Generate values from two random Gaussian vectors with different standard deviations
#' n      <- 10000 # sample size
#' sdx    <- 3 # standard deviation of the Gaussian
#' sdy    <- 5 # standard deviation of the Gaussian
#' x      <- rnorm(n, 0, sdx)
#' y      <- rnorm(n, 0, sdy)
#' # theoretical results
#' exth   <- 1 / 2 * log2((2 * pi * exp(1)) * sdx^2) # theoretical differential entropy
#' eyth   <- 1 / 2 * log2((2 * pi * exp(1)) * sdy^2) # theoretical differential entropy
#' mith   <- 0 # the theoretical mutual information is zero as x and y are independently generated
#' # Gaussian estimates
#' exg    <- entg(x)   # differential entropy for x
#' eyg    <- entg(y)   # differential entropy for y
#' mixyg  <- mig(x, y) # mutual information between x and y
#' # Kozachenko-Leonenko estimates
#' exkl   <- entkl(x)   # differential entropy for x
#' eykl   <- entkl(y)   # differential entropy for y
#' mixykl <- mikl(x, y) # mutual information between x and y
#' @export
entg <- function(x) {
  # convert x to matrix if not already
  x <- as.matrix(x)
  # get dimensions
  n <- nrow(x)
  d <- ncol(x)
  if(qr(x)$rank != d) stop("covariance matrix of x does not have full rank")
  return(1 / 2 * log2((2 * pi * exp(1))^d * det(cov(x))))
}

#' @rdname mikl
#' @export
mig <- function(x, y)
  # computation of mutual information is based on its decomposition in
  # marginal and joint differential entropies
  return(entg(x) + entg(y) - entg(cbind(x, y)))

#' @rdname mikl
#' @export
entkl <- function(x, type = "klo", k = 1) {
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

#' @rdname mikl
#' @export
mikl <- function(x, y, type = "klo", k = 1)
  # computation of mutual information is based on its decomposition in
  # marginal and joint differential entropies
  return(entkl(x, type, k) + entkl(y, type, k) - entkl(cbind(x, y), type, k))