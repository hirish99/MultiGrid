function [B] = apply_kronsum(Ax, Ay, Az, Ix, Iy, Iz, X)
    % compute b = Ax
    % A = Az x Iy x Ix + Iz x Ay x Ix + Iz x Iy x Ax
    % Ix, Iy, Iz not necessarily being identify.
    B = apply_kron(Ax, Iy, Iz, X);
    B = B + apply_kron(Ix, Ay, Iz, X);
    B = B + apply_kron(Ix, Iy, Az, X);

end;