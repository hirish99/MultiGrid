function [errors] = ksp_mg_aniso_plot(n, n_iter, ue, x_width, y_width, z_width);
    % In the anisotropic case n represents the same thing
    % there are n^3 cells
    hx=x_width/(n+1); 
    hy=y_width/(n+1);
    hz=z_width/(n+1);

    % h2i=1./(h*h); No need for h2i

    x = hx*[1:n]';
    y = hy*[1:n]';
    z = hz*[1:n]';

    e = ones(n^3,1);
    Ax = (1/hx^2)*spdiags([-e 2*e -e], -1:1, n, n);
    Ay = (1/hy^2)*spdiags([-e 2*e -e], -1:1, n, n);
    Az = (1/hz^2)*spdiags([-e 2*e -e], -1:1, n, n);
    Id = eye(n);

    A_3d = kron(Id, kron(Id, Ax)) + kron(Id, kron(Ay, Id)) + kron(Az, kron(Id, Id));

    k = [1:n]';
    Vx = sqrt(2*hx)*sin((hx*pi)*(k*k'));
    Vy = sqrt(2*hy)*sin((hy*pi)*(k*k'));
    Vz = sqrt(2*hz)*sin((hz*pi)*(k*k'));
    V_3d = kron(Vz, kron(Vy, Vx)); % Is this order correct?

    % Lam = (2*h2i)*(1-cos(h*pi*k));
    lmax = 2;
    lmin = 0.6;

    ue = V_3d*ue;

    b = A_3d*ue;
    u = 0*b;

    % M = diag(A_3d);
    % Mi = diag(1./M);

    % PCG 
    r=b;
    % z = Mi * r;
    z = vcycle_aniso(r*0, r, A_3d, n, x_width, y_width, z_width);
    p = z;
    w = A_3d*p;
    cnt = 0;
    while cnt < n_iter;
        alpha = (r' * z) / (p' * w);
        u = u + alpha * p; % Solution update
        r_old = r;
        z_old = z;
        r = r_old - alpha * (w);
        % z = Mi * r;
        z = vcycle_aniso(z*0, r, A_3d, n, x_width, y_width, z_width);
        B = (r' * z) / (r_old' * z_old);
        p = z + B * p;
        w = A_3d*p;
        cnt = cnt+1;
        errors(cnt) = norm(ue - u, Inf);

        

    end;