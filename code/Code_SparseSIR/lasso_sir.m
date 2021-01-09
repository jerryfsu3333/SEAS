function beta=lasso_sir(x,y,d,SXX)

beta= zeros(size(x,2),d);
for i=1:d
    [B,FitInfo] = lasso(x,y(:,i),'CV',5);
    ind=FitInfo.IndexMinMSE;
    beta(:,i)=B(:,ind);
end

beta=beta*(beta'*SXX*beta)^(-0.5);