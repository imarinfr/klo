from math import gamma
from numpy import mean, sqrt, cov, pi, exp, log, log2, concatenate, zeros, trace, newaxis
from numpy.linalg import matrix_rank
from scipy.linalg import sqrtm, expm, inv, det
from scipy.special import digamma
from sklearn.neighbors import KDTree


# Gaussian estimate of the differential entropy of the multivariate
# random variable x. If x is a multivariate Gaussian, then the estimate is
# asymptotically unbiased, otherwise it is just an approximation.
#
# x is a n-by-d numeric matrix, in which the n rows correspond to
# observations and the d columns to variables (or coordinates) of the
# multivariate distributions
def entg(x):
    nr, nc = x.shape
    if matrix_rank(x) != nc:
        raise Exception("covariance matrix of x does not have full rank")
    return 1 / 2 * log2((2 * pi * exp(1)) ** nc * det(cov(x.T)))


# Gaussian estimate of mutual information between the multivariate
# random variables x and y. If both x and y are multivariate Gaussians,
# then the estimate is asymptotically unbiased, otherwise it is just an
# approximation
#
# x, y are each n-by-d numeric matrices, in which the n rows correspond to
# observations and the d columns to variables (or coordinates) of the
# multivariate distributions
def mig(x, y):
    # computation of mutual information is based on its decomposition in
    # marginal and joint differential entropies
    return entg(x) + entg(y) - entg(concatenate((x, y), axis=1))


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
def entkl(x, type="klo", k=1):
    # check input
    if type not in ("kl", "klo"):
        raise Exception("incorrect type of estimator")
    if k <= 0:
        raise Exception("k-nearest neighbour 'k' must be larger than 0")
    # get dimensions
    n, d = x.shape
    if matrix_rank(x) != d:
        raise Exception("covariance matrix of x does not have full rank")
    # if type is 'kl' then calculate differential entropy with x directly.
    # If type is 'klo', then calculate differential entropy as if x were
    # Gaussian and prepare the data to calculate the offset differential
    # entropy using nearest-neighbor algorithms
    entgp = 0
    if type == "klo":
        entgp = entg(x)
        x = x.dot(inv(sqrtm(cov(x.T)))) / sqrt(2 * pi * exp(1))
    # the nearest-neighbour distances
    knn = KDTree(x, metric="euclidean")
    dist, idx = knn.query(x, k=k + 1)
    dist = dist[:, k]
    ldist = zeros(n)  # init log distances with zeros
    ldist[dist != 0] = log(dist[dist != 0])  # log 0 = 0
    # differential entropy calculation in nats
    ent = d * mean(ldist) + log(2 * pi ** (d / 2) / (d * gamma(d / 2))) + log(n - 1) - digamma(k)
    # return differential entropy in bits
    return log2(exp(ent)) + entgp


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
def mikl(x, y, type="klo", k=1):
    # computation of mutual information is based on its decomposition in
    # marginal and joint differential entropies
    return entkl(x, type, k) + entkl(y, type, k) - entkl(concatenate((x, y), axis=1), type, k)