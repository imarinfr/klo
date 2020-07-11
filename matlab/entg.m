function ent = entg(x)
% ENTG: Gaussian estimate of the differential entropy of the multivariate
% random variable x. If x is a multivariate Gaussian, then the estimate is
% asymptotically unbiased, otherwise it is just an approximation.
%
% x is a n-by-d numeric matrix, in which the n rows correspond to
% observations and the d columns to variables (or coordinates) of the
% multivariate distributions

% check input
if nargin < 1, error('please revise input'); end

% get dimension
[~, d] = size(x);
if rank(x) ~= d
    error('covariance matrix of x does not have full rank');
end

ent = 1 / 2 * log2((2 * pi * exp(1))^d * det(cov(x)));