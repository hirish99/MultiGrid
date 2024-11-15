function [u] = vcycle_aniso(u, rhs, A, n, x_width, y_width, z_width);
    u = jac_smooth(u, rhs, A);
    % u = cheb_smooth(u, rhs, A);

    r = rhs - A*u;

    % restriction
    nc = (n+1)/2 - 1;

    x = [1:n]'/(n+1) * x_width;
    xc = [1:nc]'/(nc+1) * x_width;

    y = [1:n]'/(n+1) * y_width;
    yc = [1:nc]'/(nc+1) * y_width;

    z = [1:n]'/(n+1) * z_width;
    zc = [1:nc]'/(nc+1) * z_width;

    J_1dx = lin_interp_mat([0; x; 1],[0; xc; 1]);
    J_1dy = lin_interp_mat([0; y; 1],[0; yc; 1]);
    J_1dz = lin_interp_mat([0; z; 1],[0; zc; 1]);

    J_1dx = J_1dx(2:end-1,2:end-1);
    J_1dy = J_1dy(2:end-1,2:end-1);
    J_1dz = J_1dz(2:end-1,2:end-1);

    J = kron(J_1dz, kron(J_1dy, J_1dx)); %Is this order correct?

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