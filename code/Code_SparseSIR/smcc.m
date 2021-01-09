function  s = smcc(A, B, X)
% smcc: Square Multiple correlation coefficient
SXX = cov(X);
a = A'*SXX*A;
a = invmat(a);
b = A'*SXX*B;
c = B'*SXX*B;
c = invmat(c);
d = B'*SXX*A;
m = a*b*c*d;
s = trace(m);
end