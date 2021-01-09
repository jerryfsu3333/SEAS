function zz = tab(Mod,n,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function makes talbe for different models with various (n,p) setting.
% Input: Mod specify the model. n,p are the sample size and dimension. 
% Output: a 20* 16 matrix, whose odd columns is just like Table 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2;               % number of repeating times
s = floor(sqrt(p));   % number of active variables
rho = 0.1;            % penalty parameter for ADMM algorithm
ind = randi(p, 1, s); % the index of the active variables

% Models 1-4 are single index model, Model 5's structural dimension is 2.
if Mod<5
    d = 1;
else
    d = 2;
end

% Generate coefficients for the active variables
b = zeros(p, d);
b(ind,:) = randi([-2 2],s,d);
h = 10;               % number of slices
eta = 2;
e = zeros(k, 3, 4);   
r = e;
e2 = zeros(k, 4);     % general loss
r2= e2;               % correlation loss
%% Generate 4 different covariance matrices
sig = zeros(p, p, 4);
sig(:,:,1) = eye(p);
sig(:,:,2) = 0.6*eye(p)+0.4*ones(p,p);
sig(:,:,3) = toeplitz(0.5.^(0:1:(p-1)));
a = eye(p);
for i = 1:p
    for j = 1:p
        if abs(i-j)==1
            a(i,j) = 0.5;
        elseif abs(i-j)==2
            a(i,j) = 0.4;
        end
    end
end
ts = inv(a);
aa = zeros(p);
for i=1:p
    for j=1:p
        aa(i,j) = ts(i,j)/sqrt(ts(i,i)*ts(j,j));
    end
end
sig(:,:,4) = aa;
%% Generate Models 1 - 5
bx = zeros(p,d);
pb = zeros(p,p);
for l = 1:4             % 4 types of covariance matrices
    covx = sig(:,:, l);
    t1 = b(ind,:);
    bx(ind,:) = t1/sqrtm(t1'*covx(ind,ind)*t1);
    bx = real(bx);
    pb(ind,ind) = t1/(t1'*covx(ind,ind)*t1)*t1';
    for i = 1:k         % Repeat the simulation k times
        x=mvnrnd(zeros(p,1),covx,n);
        SXX = cov(x);
        switch Mod
            case 1
                y = x*b + sin(x*b) + randn(n,1);
            case 2
                y = 2*atan(x*b) + randn(n,1);
            case 3
                y = (x*b).^3 + randn(n,1);
            case 4
                y = asinh(x*b) + randn(n,1);
            case 5
                y=exp(x*b(:,1)).*sign(x*b(:,2)) + 0.2*normrnd(0,1,n,1);
        end
        [M,Exy,ynew]= fun_sir(x, y, h, d);
		% DT-SIR estimator 
        b1 = dtsir(x, y, h, d);
		% Natural estimator 
        b2 = ssir_natural(SXX, M, n, d, rho, eta, 1e-4, 1e2);
		% Refined estimator
        b3 = ssir_refine(x, Exy, SXX, b2, p, d);
		
        e(i, 1, l) = norm(b1*b1'-pb,'fro')^2;
        e(i, 2, l) = norm(b2*b2'-pb,'fro')^2;
        e(i, 3, l) = norm(b3*b3'-pb,'fro')^2;
        r(i, 1, l) = 1-smcc(b1,b,x)/d;
        r(i, 2, l) = 1-smcc(b2,b,x)/d;
        r(i, 3, l) = 1-smcc(b3,b,x)/d;
        
		% Lasso-SIR estimator
        b4=lasso_sir(x,ynew,d,SXX);
        if any(isnan(b4(:))) ==1
            e2(i, l) = 100;
            r2(i, l) = 100;
        else
            e2(i, l) = norm(b4*b4'-pb,'fro')^2;
            r2(i, l) = 1-smcc(b4,b,x)/d;
        end
    end
end
%% Collect the simulation results
aa = zeros(4,2);
bb = aa;
for l = 1:4
    t1= e2(:,l);
    t2= r2(:,l);
    t1(t1==100)=[];
    t2(t2==100)=[];
    aa(l,:) = [mean(t1), std(t1)];
    bb(l,:) = [mean(t2), std(t2)];
end
%%
z = zeros(4, 12);
for l = 1: 4
    z(l, 1:2:12) = [mean(e(:,:, l), 1), mean(r(:,:, l), 1)];
    z(l, 2:2:12) = [std(e(:,:, l), 1), std(r(:,:, l), 1)];
end
zz = [z(:,1:2), aa, z(:,3:8), bb, z(:,9:12)];
%name = strcat('save_tab_M',num2str(Mod),'_',num2str(n),'_',num2str(p));
%save(name)