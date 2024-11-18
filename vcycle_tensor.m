function [U] = vcycle_tensor(U, Rhs, Ax, Ay, Az, Ix, Iy, Iz, nx, ny, nz, Lx, Ly, Lz);

    AU = apply_kronsum(Ax, Ay, Az, Ix, Iy, Iz, U);
    R = Rhs - AU;
    
    % jacobi smoothing
    nsmooth = 4;
    omega = 2/3;
    d = (Ax(1,1).*Iy(1,1).*Iz(1,1) + Ix(1,1).*Ay(1,1).*Iz(1,1) + Ix(1,1).*Iy(1,1).*Az(1,1));
    di = 1./d;

    for nu=1:nsmooth;
        U=U+omega.*di.*R;
        AU = apply_kronsum(Ax, Ay, Az, Ix, Iy, Iz, U);
        R = Rhs - AU;
    end;

    % restriction
    ncx = (nx+1)/2 - 1;
    ncy = (ny+1)/2 - 1;
    ncz = (nz+1)/2 - 1;

    x = [1:nx]'/(nx+1) * Lx;
    xc = [1:ncx]'/(ncx+1) * Lx;
    y = [1:ny]'/(ny+1) * Ly;
    yc = [1:ncy]'/(ncy+1) * Ly;
    z = [1:nz]'/(nz+1) * Lz;
    zc = [1:ncz]'/(ncz+1) * Lz;

    Jx = lin_interp_mat([0; x; Lx],[0; xc; Lx]);
    Jy = lin_interp_mat([0; y; Ly],[0; yc; Ly]);
    Jz = lin_interp_mat([0; z; Lz],[0; zc; Lz]);

    Jx = Jx(2:end-1,2:end-1);
    Jy = Jy(2:end-1,2:end-1);
    Jz = Jz(2:end-1,2:end-1);

    Acx = Jx' * Ax * Jx;
    Acy = Jy' * Ay * Jy;
    Acz = Jz' * Az * Jz;
    
    JJx = Jx' * Ix * Jx;
    JJy = Jy' * Iy * Jy;
    JJz = Jz' * Iz * Jz;


    if ncx == 1;
        Jx = full(Jx);
        Jy = full(Jy);
        Jz = full(Jz);
    end;
    
    Rc = apply_kron(Jx',Jy',Jz',R);

    % We are solving A x= rc
    if ncz == 1 || ncy == 1 || ncx == 1;
        Ec = (JJz .* JJy .* Acx + JJz .* Acy .* JJx + Acz .* JJy .* JJx) \ Rc;
    else
        Ec = vcycle_tensor(Rc*0, Rc, Acx, Acy, Acz, JJx, JJy, JJz, ncx, ncy, ncz, Lx, Ly, Lz);
    end;
    % coarse grid correction
    E = apply_kron(Jx, Jy, Jz, Ec);

    U = U + E;
    AU = apply_kronsum(Ax, Ay, Az, Ix, Iy, Iz, U);
    R = Rhs - AU;

end;