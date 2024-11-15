clear;
hdr;

k = 5;
% In the anisotropic case n represents the same thing
% there are n^3 cells
x_width = 2;
y_width = 1;
z_width = 1;
n=8; 
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
V_3d = kron(Vz, kron(Vy, Vx)); %Is this order correct?

k = [1:n]';
V = sqrt(2*h)*sin((h*pi)*(k*k'));
V_3d = kron(V, kron(V, V));

Lam = (2*h2i)*(1-cos(h*pi*k));
lmax = 2;
lmin = 0.6;

ue = rand(n^3,1);
ue = V_3d*ue;

b = A_3d*ue;
u = 0*b;

% M = diag(A_3d);
% Mi = diag(1./M);

% PCG
r=b;
% z = Mi * r;
z = vcycle(r*0, r, A_3d, n);
p = z;
w = A_3d*p;
cnt = 0;
while norm(r) > 1e-8;
    alpha = (r' * z) / (p' * w);
    u = u + alpha * p; % Solution update
    r_old = r;
    z_old = z;
    r = r_old - alpha * (w);
    % z = Mi * r;
    z = vcycle(z*0, r, A_3d, n);
    B = (r' * z) / (r_old' * z_old);
    p = z + B * p;
    w = A_3d*p;

    er = ue - u;
    norm(er)
    cnt = cnt+1;

end;

cnt