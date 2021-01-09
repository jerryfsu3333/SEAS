function b2 = ssir_refine(x, Exy, SXX, b1, p, d)

% Inputs:
% =======
% x, Exy, SXX: X; Cov[E(X|Y)]; Cov(X);
% b1, d: initial estimators for b; structural dimension
% 
% Outputs:
% ========
% b2: estimators for b
%ssir_refine_cv(x, Exy, SXX, b(:,:,2), p, d);
%b1 = b(:,:,3);

% n = size(x,1);
% x = x - ones(n,1) * mean(x);
% Exy = Exy - ones(n,1) * mean(Exy);

Y = Exy*b1;
B = zeros(p+1,d);
if (d==1)
    fit = glmnet(x, Y);
    cvfit = cvglmnet(x, Y);
    B=glmnetCoef(fit,cvfit.lambda_min,true);
else 
    fit = glmnet(x, Y, 'mgaussian');
    cvfit = cvglmnet(x, Y, 'mgaussian');
    ncoef=glmnetCoef(fit,cvfit.lambda_min,true);
    for i = 1:d
        B(:,i) = ncoef{i};
    end
    %B = cell2mat(ncoef);
end
B(1,:) = [];
ind = sum(abs(B),2)~=0;
b2 = zeros(p, d);
b2(ind,:) = B(ind,:)/sqrtm(B(ind,:)'*SXX(ind, ind)*B(ind,:));
% b2 = b2/sqrtm(b2'*SXX*b2);