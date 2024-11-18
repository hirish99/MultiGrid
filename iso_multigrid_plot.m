function [errors] = iso_multigrid_plot(n, n_iter);

    h=1/(n+1); h2i=1./(h*h);

    x = h*[1:n]';
    e = ones(n^3,1);
    A = spdiags([-e 2*e -e], -1:1, n, n);

    A_2d = kron(A, speye(n)) + kron(speye(n), A);
    A_3d = kron(A, speye(n^2)) + kron(speye(n), A_2d);

    A_3d = h2i * A_3d;

    k = [1:n]';
    V = sqrt(2*h)*sin((h*pi)*(k*k'));
    V_3d = kron(V, kron(V, V));

    Lam = (2*h2i)*(1-cos(h*pi*k));
    lmax = 2;
    lmin = 0.6;




    % ue=0*x;
    ue=rand(n^3,1);
    ue=V_3d*e;

    b = A_3d*ue;
    u=0*b;

    r=b;

    cnt = 0;
    r_norm = 1;

    while cnt < n_iter;
        u = vcycle(u, b, A_3d, n);
        % residual = b - A_3d*u;
        cnt = cnt+1;
        errors(cnt) = norm(u - ue, Inf);
    end;