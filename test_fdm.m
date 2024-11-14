m = 256;
L = 1
F = rand(m,m,m);

h = 1/(m+1);

e = ones(m,1);
A = spdiags([-e 2*e -e], -1:1, m, m) ./ h^2;
X = fdm(m, L, F);

AX = tensorprod(X,A,3,1);
AX = AX + permute(tensorprod(X,A,2,1), [1,3,2]);
AX = AX + tensorprod(A,X,1,1);


norm(reshape(AX-F, m^3,1))