function sol = ssir_natural(SXX, M, n, d, rho1, eta, tol, maxiter)

% Inputs:
% =======
% SXX, M:   Covariance matrices Cov(X), Cov[E(X|Y)];
% d:        nuclear norm constraint
% rho1:     penalty parameter on the l_1 norm of the solution, scaled by sqrt(log(max(p1,p2))/n)
% eta = 2;  parameter in the augmented lagrangian
% tol:      tolerance level for convergence in ADMM
% maxiter:  maximum number of iterations in ADMM
%
% Outputs:
% ========
% sol:     optimum of the convex program
% sol1:    resacled optimum of the convex program, premultiplied by
p = size(SXX, 1);
%% compute the root inverse of SXX
[VX, DX] = eig(SXX);
DX = diag(DX);
idx = (abs(DX) > max(abs(DX)) * 1e-6); %find the index that eigenvalue not equals to 0
DX = sqrt(DX(idx));
SXroot = VX(:,idx) * diag(DX) * (VX(:,idx)');
SXrootInv = VX(:,idx) * diag(1./DX) * (VX(:,idx)');

B = SXrootInv * M * SXrootInv;

A = @(varargin)linop_crosprod(p, p, SXroot, SXroot, varargin{:} );

x_cur = 0;                        % G0 in Algorithm 1 of the paper.
y_cur = zeros(p);             % H0 in Algorithm 1 of the paper.

%[U, D, V] = svd(M, 'econ');
[V, D] = eigs(M);
d1 = diag(D);
t = cap_soft_th(d1, d, tol);
z_cur = V * diag(t) * V';         % F0 in Algorithm 1 of the paper.

niter = 0;
initer = 0;

opts = [];
opts.printEvery = Inf;
opts.maxIts = 25;

while 1
    niter = niter + 1;
    initer = initer + 1;
    % Update F
    z_old = z_cur;
    Temp = x_cur - y_cur ./ eta +  B ./ eta;
    z_cur = tfocs(smooth_quad, {A, -Temp}, prox_l1(2 * rho1 * sqrt(log(p)/n) / eta), z_old, opts); % update F
    
    % Update G
    x_old = x_cur;
    Temp = y_cur ./ eta + A(z_cur, 1);  % A(z_cur, 1) is SXroot F SYroot
    %[U, D, V] = svd(Temp, 'econ');
    [V, D] = eigs((Temp + Temp')/2);
    d2 = diag(D);
    t = cap_soft_th(d2, d, tol);
    x_cur = V * diag(t) * V';
    
    % Update H in eq(117)
    y_old = y_cur;
    y_cur = y_old + eta .* (A(z_cur,1) - x_cur); % H update in Algorithm 1 of the paper.
    
    if max(eta*sqrt(log(p)/n)*norm(x_cur - x_old, 'fro'), norm(z_cur - z_old, 'fro')) < tol
        break
    end
    if niter == maxiter
        %disp('Maximum number of iterations reached.');
        break
    end
end

[vdr,~]=eigs((z_cur+z_cur')/2);
vdr=real(vdr);
sol = vdr(:, 1:d);
sol = sol*(sol'*SXX*sol)^(-0.5);