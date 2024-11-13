function [u] = cheb_smooth(u, rhs, A);
    r = rhs - A*u;
    nsmooth = 4;
    omega = 2/3;
    di = 1./A(1,1);

    % TODO: add cheb_smooth

end;