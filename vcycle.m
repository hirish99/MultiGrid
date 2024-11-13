function [u] = vcycle(u, rhs, A, n);
    u = jac_smooth(u, rhs, A);
    % u = cheb_smooth(u, rhs, A);

    r = rhs - A*u;

    % restriction
    nc = (n+1)/2 - 1;

    x = [1:n]'/(n+1);
    xc = [1:nc]'/(nc+1);
    J_1d = lin_interp_mat([0; x; 1],[0; xc; 1]);
    J_1d = J_1d(2:end-1,2:end-1);
    J = kron(J_1d, kron(J_1d, J_1d));

    rc = J' * r;

    ec = zeros(size(rc));
    Ac = J'*A*J;
    % We are solving A x= rc
    if size(rc) == 1;
        ec = Ac \ rc;
    else
        nc = (n+1)/2-1;
        ec = vcycle(ec, rc, Ac, nc);
    end;
    % coarse grid correction
    ef = J * ec;
    u = u + ef;
    r = rhs - A*u;

    % e=ue-u;


end;