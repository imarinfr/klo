function ent = entkl(x, type, k)
% ENTKL: The Kozachenko-Leonenko (KL) estimate of the differential entropy
% of the multivariate random variable x.
%
% If type is 'kl', the raw KL estimator is computed.
% If type is 'klo', the offset KL estimator is computed.
%
% k is the nearest neighbor for which to search, 1 the nearest neighbor, 2,
% the second nearest, and so on
%
% x is a n-by-d numeric matrix, in which the n rows correspond to
% observations and the d columns to variables (or coordinates) of the
% multivariate distributions

% check input
if nargin < 1, error('please revise input'); end
if nargin < 2, type = 'klo'; end
if nargin < 3, k = 1; end
if ~strcmp(type, 'kl') && ~strcmp(type, 'klo')
    error('incorrect type of estimator');
end
assert(k > 0, 'k-nearest neighbour ''k'' must be larger than 0')
% get dimensions
[n, d] = size(x);
if rank(x) ~= d
    error('covariance matrix of x does not have full rank');
end
% if type is 'kl' then calculate differential entropy with x directly.
% If type is 'klo', then calculate differential entropy as if x were
% Gaussian and prepare the data to calculate the offset differential
% entropy using nearest-neighbor algorithms
entgp = 0;
if strcmp(type, 'klo')
    entgp = entg(x);
    x     = x / sqrtm(cov(x)) / sqrt(2 * pi * exp(1));
end
% the nearest-neighbour search
[~, dist] = knnsearch(x, x, 'K', k + 1);
ldist = log(dist(:,k + 1));
ldist(ldist == -Inf) = 0; % log 0 = 0
% differential entropy calculation in nats
ent = d * mean(ldist) + log(2 * pi^(d / 2) / (d * gamma(d / 2))) + log(n - 1) - psi(k);
% return differential entropy in bits
ent = log2(exp(ent)) + entgp;