function Be = dtsir(x, y, h, d)
% This is the main function that achives DT_SIR algorithm
%    USAGE:
%  - outputs:
%    Be: the estimate of Bx, which satisfy Bx'* Sigma_x * Bx = I_d;

%  - inputs:
%    y: response vector.
%    x: predictors matrix.
%    d: the dimension of the central subspace.
%    h: number of slice.
%% threshold
[n, p] = size(x);
SXX = cov(x);
%repeat k times
k = 5;
Z = cell(1,k);
thres = zeros(1, k);
for i = 1:k
    Z{i} = mvnrnd(zeros(p,1), eye(p), n);
    [SM_Z, ] = sir(Z{i}, y, h);
    thres(i) = max(diag(SM_Z));
end
threshold = mean(thres);
%% screening
[SM, ~] = sir(x, y, h);
varH = diag(SM);
active = find(varH >threshold);
l = length(active);
if isempty(active)
    Be= zeros(p, d);
else
    xs = x(:, active);              % x after selection
    zs = zscore(xs);
    SMs = sir(zs, y, h);          % Sample kernal matrix after selection
    [Vs,~] = eigs(SMs);
    if l<d
        Bs = [Vs, ones(l, d-l)];
        Be = zeros(p,d);
        Be(active, :) = Bs;
    else
        Bs = Vs(:, 1:d);
        Be = zeros(p,d);
        Be(active, :) = Bs;
        Be = Be/sqrtm(Be'* SXX * Be);
    end
end