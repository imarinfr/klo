function mi = mig(x, y)
% MIG: Gaussian estimate of mutual information between the multivariate
% random variables x and y. If both x and y are multivariate Gaussians,
% then the estimate is asymptotically unbiased, otherwise it is just an
% approximation
%
% x, y are each n-by-d numeric matrices, in which the n rows correspond to
% observations and the d columns to variables (or coordinates) of the
% multivariate distributions

% check input
if nargin < 2, error('please revise input'); end

% computation of mutual information is based on its decomposition in
% marginal and joint differential entropies
mi = entg(x) + entg(y) - entg([x y]);