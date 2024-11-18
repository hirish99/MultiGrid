function [B] = apply_kron(Ax, Ay, Az, X)
    % compute vec(B) = (Az kron Ay kron Ax) * vec(X)
    B = tensorprod(Ax, X, 2, 1);
    B = permute(tensorprod(Ay, B, 2, 2), [2,1,3]);
    B = tensorprod(B, Az, 3, 2);
end;