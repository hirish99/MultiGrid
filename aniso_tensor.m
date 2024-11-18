clear;
hdr;

Lx = 1;
Ly = 1;
Lz = 1;

% nx, ny, nz must be equal for now.
nx=31;
ny=31;
nz=31;

hx=Lx/(nx+1);
hy=Ly/(ny+1);
hz=Lz/(nz+1);

x = hx*[1:nx]';
y = hy*[1:ny]';
z = hz*[1:nz]';
e = ones(nx^2,1);
Ax = spdiags([-e 2*e -e], -1:1, nx, nx) / (hx^2);
Ay = spdiags([-e 2*e -e], -1:1, ny, ny) / (hy^2);
Az = spdiags([-e 2*e -e], -1:1, nz, nz) / (hz^2);


kx = [1:nx]';
ky = [1:ny]';
kz = [1:nz]';
Vx = sqrt(2*hx)*sin((hx*pi)*(kx*kx'));
Vy = sqrt(2*hy)*sin((hy*pi)*(ky*ky'));
Vz = sqrt(2*hz)*sin((hz*pi)*(kz*kz'));


Ue=rand(nx,ny,nz);
% ue = reshape(Ue, [], 1);
Ue = apply_kron(Vx, Vy, Vz, Ue);

% A_2d = kron(A, speye(n)) + kron(speye(n), A);
% A_3d = kron(A, speye(n^2)) + kron(speye(n), A_2d);
% V_3d = kron(V, kron(V, V));
% ue = V_3d*ue;

Ix = eye(nx);
Iy = eye(ny);
Iz = eye(nz);
B = apply_kronsum(Ax,Ay,Az,Ix,Iy,Iz,Ue);


cnt = 0;
U=0*B;
R_norm = 1;
while R_norm > 1e-8;
% for kk = 1:10;
    U = vcycle_tensor(U, B, Ax, Ay, Az, Ix, Iy, Iz, nx, ny, nz, Lx, Ly, Lz);
    AU = apply_kronsum(Ax, Ay, Az, Ix, Iy, Iz, U);
    R = B - AU;
    cnt = cnt+1;
    R_norm = norm(reshape(full(R), [] ,1));
    disp(R_norm);
end;


disp(norm(reshape(full(U - Ue), [], 1)))
disp(cnt)
