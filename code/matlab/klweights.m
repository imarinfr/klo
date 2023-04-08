% This is a lightly edited copy of Berrett's L2OptW function in IndepTest
% R package. It calculates a weight vector to be used for the weighted
% Kozachenko–Leonenko estimator. The weight vector has minimum L_2 norm
% subject to the linear and sum-to-one constraints of (2) in Berrett,
% Samworth and Yuan, Annals of Statistics (2019).
function w = klweights(k, d)
dprime = floor(d / 4);
if k == 1
    w = 1;
elseif dprime == 0
    w = [zeros(1, k - 1), 1];
else
    g = zeros(dprime + 1, k);
    g(1,:) = ones(1, k);
    for l = 1:dprime
        g(l + 1,:) = exp(gammaln((1:k) + 2 * l / d) - gammaln(1:k));
    end
    w = [1, zeros(1, dprime)] / (g * g') * g;
end