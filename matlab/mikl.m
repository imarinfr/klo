function mi = mikl(x, y, type, k)
% MIKL: The Kozachenko-Leonenko (KL) estimate of mutual information between
% the multivariate random variables x and y. 
%
% If type is 'kl', the raw KL estimator is computed.
% If type is 'klo', the offset KL estimator is computed.
%
% k is the nearest neighbor for which to search, 1 the nearest neighbor, 2,
% the second nearest, and so on
%
% x, y are each n-by-d numeric matrices, in which the n rows correspond to
% observations and the d columns to variables (or coordinates) of the
% multivariate distributions

% check input
if nargin < 2, error('please revise input'); end
if nargin < 3, type = 'klo'; end
if nargin < 4, k = 1; end
% computation of mutual information is based on its decomposition in
% marginal and joint differential entropies
mi = entkl(x, type, k) + entkl(y, type, k) - entkl([x y], type, k);