function [X] = fdm(m, L, F);
    % m: number of grid points in one dimension
    % L: domain length
    % F: rhs, shape is (m,m,m)
    h = L/(m+1);
    i= [1:m]';
    ij = i*i';
    scale = sqrt(2*h);
    S = scale * sin(ij * (pi/(m+1)));
    eigvals = 2 ./ h^2 * (1 - cos(pi .* i / (m+1)));
    eigvals_3d = (kron(ones(m*m,1), eigvals) + kron(kron(ones(m,1), eigvals), ones(m,1)) + kron(eigvals, ones(m*m,1)));
    eigvals_3d = reshape(eigvals_3d,m,m,m);



    F_hat = tensorprod(F, S', 3 ,1);
    F_hat = permute(tensorprod(F_hat, S', 2,1), [1,3,2]);
    F_hat = tensorprod(S', F_hat, 1,1);
    X_hat = F_hat ./ eigvals_3d;
    X = tensorprod(X_hat, S, 3,1);
    X = permute(tensorprod(X, S, 2, 1), [1,3,2]);
    X = tensorprod(S, X, 1, 1);

end;

