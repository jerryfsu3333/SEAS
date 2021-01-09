%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This m-file produces Fig 1, 2(in the paper), Table 1, 2, 3(in Appendix)
% z is a 20*16 matrix, whose odd columns are the mean of the errors, and the even columns are the standard error. 
% The submitted paper only presents the mean of errors due to space limitation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc
tic;
warning('off');
%%
p = 500;                 % dimension of variables
Mod = 1;                 % Model 
k = 2;                   % number of repeating times
t = 10;                  % length of sequence for scaled sample size
s = floor(sqrt(p));      % number of active variables
ind = randi(p, 1, s);    % the index of the active variables

% Models 1-4 are single index model, Model 5's structural dimension is 2.
if Mod<5
    d = 1;
else
    d = 2;
end

% Generate coefficients for the active variables
b = zeros(p, d);
b(ind,:) = randi([-2 2],s,d); 
eta = 2;
rho = 0.2;                    % penalty parameter for ADMM algorithm
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
%% Running our two step algorithm
[e1, e2, r1, r2]= deal(zeros(t, k, 4));
bx = zeros(p,d);
pb = zeros(p,p);
for l = 1: 4
    covx = sig(:,:,l);
    t1 = b(ind,:);
    bx(ind,:) = t1/sqrtm(t1'*covx(ind,ind)*t1);
    bx = real(bx);
    pb(ind,ind) = t1/(t1'*covx(ind,ind)*t1)*t1';
    for j = 1:k
        for i = 1:t
            kap = 3*i;
            n = floor(kap*(s*log(exp(1)*p/s)));
            if n > 1000
                h = 20;
            else
                h = 10;
            end
            x = mvnrnd(zeros(p,1), covx, n);
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
            [M, Exy] = sir(x, y, h);
            
            % natural estimator 
            b1 = ssir_natural(SXX, M, n, d, rho(1), eta, 1e-4, 2e2);
            % refined estimator
            b2 = ssir_refine(x, Exy, SXX, b1, p, d);
            
            e1(i, j, l) = norm(b1*b1'-pb,'fro')^2;
            r1(i, j, l) = 1-smcc(b1,b,x)/d;
            e2(i, j, l) = norm(b2*b2'-pb,'fro')^2;
            r2(i, j, l) = 1-smcc(b2,b,x)/d;
        end
    end
end
%% Collect Results for plotting
z1 = zeros(t, 2, 4);
z2 = z1;
for l = 1:4
    z1(:,:, l) = [mean(e1(:,:,l), 2), mean(e2(:,:,l), 2)];
    z2(:,:, l) = [mean(r1(:,:,l), 2), mean(r2(:,:,l), 2)];
end
ch = {'-.^b', '--or'};
kappa = 3*(1:t);
time = toc
%% Plot of General Loss
figure
tit{1} = '$\mathbf{\Sigma}_1= \mathbf{I}_p$';
tit{2} = '$\mathbf{\Sigma}_2: \sigma_{ij} =1(i=j) + 0.6 \cdot 1(i\neq j)$';
tit{3} = '$\mathbf{\Sigma}_3: \sigma_{ij} = 0.5^{|i-j|}$';
tit{4} = '$\mathbf{\Sigma}_4: \mathrm{~SparseInv}$';
for l = 1: 4
    subplot(2, 2, l)
    for i = 1:2
        plot(kappa, z1(:, i, l), ch{i})
        hold on
    end
    legend({'Natural $\hat{\beta}^{\star}$', 'Refined $\tilde{\beta}^{\star}$'},'Interpreter','latex', 'Location', 'NorthEast')
    xlabel({'$t$'},'Interpreter','latex')
    ylabel({'$\hat{E}L_G(\hat\beta,\beta)$'},'Interpreter','latex');
    te = title(tit{l});
    set(te,'Interpreter','latex');
end
c1 = strcat('M',num2str(Mod),'_1.eps');
print(c1,'-depsc2', '-r600');
%% Plot of Correlation Loss
figure
for l = 1: 4
    subplot(2, 2, l)
    for i = 1:2
        plot(kappa, z2(:, i, l), ch{i})
        hold on
    end
    legend({'Natural $\hat{\beta}^{\star}$', 'Refined $\tilde{\beta}^{\star}$'},'Interpreter','latex', 'Location', 'NorthEast')
    xlabel({'$t$'},'Interpreter','latex')
    ylabel({'$\hat{E}L_{\rho}(\hat{\beta}, \beta)$'},'Interpreter','latex');
    te = title(tit{l});
    set(te,'Interpreter','latex');
end
c2 = strcat('M',num2str(Mod),'_2.eps');
print(c2,'-depsc2', '-r600');