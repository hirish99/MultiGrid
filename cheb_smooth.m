function [u] = cheb_smooth(u, rhs, A);
    r = rhs - A*u;


    d = diag(A);
    di=1./d;
    lmax = 5.9;
    lmin = 0.7;

    theta = .5*(lmax+lmin);
    delta = .5*(lmax-lmin);
    sigma = theta/delta;
    rho0  = 1./sigma;
    r=di.*r;
    d=r/theta;


    for iter = 1:4;
        u = u + d;
        r = r - di.*(A*d);
        rho1=rho0;
        rho1 = 1/(2*sigma-rho0);
        d = (rho1*rho0)*d + (2*rho1/delta)*r;
    end;


end;