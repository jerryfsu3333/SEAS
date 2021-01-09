function [A, B, C] = invmat(X)
% computing the inverse of matrix X, i.e. X^(-1), X^(-0.5)
% A: inverse of X
% B: root of X
% C: root inverse of X
%[n, p] = size(X);
[VX, DX] = eig(X);
DX = diag(DX);
idx = (abs(DX) > max(abs(DX)) * 1e-6); %find the index that eigenvalue not equals to 0
%idx = (abs(DX) > 2*sqrt(log(p)/n));  % Precision matrix estimation by Bickle and Levina (2008, AOS)
DX = DX(idx);
A = VX(:,idx) * diag(1./DX) * (VX(:,idx)');
B = VX(:,idx) * diag(sqrt(DX)) * (VX(:,idx)');
C = VX(:,idx) * diag(1./sqrt(DX)) * (VX(:,idx)');
end