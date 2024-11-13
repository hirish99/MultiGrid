function [u] = jac_smooth(u, rhs, A);
    r = rhs - A*u;
    nsmooth = 4;
    omega = 2/3;
    di = 1./A(1,1);

    for nu=1:nsmooth;
        u=u+omega*di*r;
        r=rhs-A*u;
    end;

end;