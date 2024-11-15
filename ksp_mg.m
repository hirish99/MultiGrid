clear;
hdr;

k = 5;
n=2^k-1; h=1/(n+1); h2i=1./(h*h);

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