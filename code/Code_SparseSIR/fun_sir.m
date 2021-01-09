function [M,Exy,ynew]= fun_sir(x, y, h, d)
%This is the main function that achives SIR algorithm
%    USAGE:
%  - outputs:
%    SM: sample kernal matrix Cov[E(x|y)];
%    Exy: sample of E(x|y);

%  - inputs:
%    y: response vector.
%    x: predictors matrix.
%    h: number of slice.
[n, p] = size(x);
x = x-ones(n,1)*mean(x);
y1 = unique(y);
nc = length(y1);
if nc < n/2
    h = nc;
else
    h = h;
end
m = zeros(h,p);
Exy = zeros(n, p);
a = zeros(n,h);
yt = zeros(n,1);
ph = zeros(1,h);
if nc < n/2
    for i = 1:h
        ind = find(y==y1(i));
        ni = length(ind);
        ph(i)=ni/n;
        yt(ind) = i;
        a(:,i) = (yt==i);
        xi = x(ind, :);
        m(i,:) = mean(xi);
        Exy(ind,:) = repmat(m(i,:), ni, 1);
    end
else
    [~,YI] = sort(y);
    c = round(n/h);
    for i = 1:h
        if (i < h)
            ind = YI( ((i-1)*c + 1):(i*c) );
        else
            ind = YI( ((h-1)*c + 1):n );
        end
        yt(ind) = i;
        ni = length(ind);
        ph(i) = ni/n;
        a(:,i) = (yt==i);
        xi = x(ind, :);
        m(i,:) = mean(xi);
        Exy(ind,:) = repmat(m(i,:), ni, 1);
    end
end
M = x'* a * diag(1./ph) * a' * x/(n^2);
M = (M + M')/2;
[V,D]=eigs(M,d);
V=real(V);
D=real(D);
ynew = 1/n*a * diag(1./ph) * a' * x*V/D;