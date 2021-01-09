function [SM,  Exy]= sir(X, Y, h)
% This is the main function that achives SIR algorithm
%    USAGE:
%  - outputs:
%    SM: sample kernal matrix Cov[E(x|y)]
%    Exy: sample of E(x|y);

%  - inputs:
%    Y: response vector.
%    X: predictors matrix.
%    h: number of slice.

[n, p] = size(X);

[~, y.index] = sort(Y);
YI = y.index;
Hn = round(n/h);
m = zeros(h,p);
ph = zeros(1, h);
xmean = mean(X);
Exy = zeros(n, p);
for i = 1:h
    if (i < h)
        ind = YI( ((i-1)*Hn + 1):(i*Hn) );
    else
        ind = YI( ((h-1)*Hn + 1):n );
    end
    ni = length(ind);
    ph(i) = ni/n;
    xi = X(ind, :);
    m(i,:) = mean(xi);
    Exy(ind,:) = repmat(m(i,:), ni, 1);
    m(i,:) = m(i,:) - xmean;
end
p = diag(ph);
SM = m' * p * m;
%SM = (SM + SM')/2;          % symmetry
end