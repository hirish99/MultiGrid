hdr;

n=31; h=1/(n+1); h2i=1./(h*h);

x = h*[1:n]';
e = ones(n^3,1);
A = spdiags([-e 2*e -e], -1:1, n, n);

A_2d = kron(A, eye(n)) + kron(eye(n), A);
A_3d = kron(A, eye(n^2)) + kron(eye(n), A_2d);

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



di=1./A_3d(1,1);

% PCG
r=b;
% z = Mi*r;
z = r*0;
z = vcycle(z, r, A_3d, n);
p = z;
w = A_3d*p;
for k=1:30;
    alpha = (r' * z) / (p' * w);
    u = u + alpha * p; % Solution update
    r_old = r;                                
    z_old = z;                                     
    r = r_old - alpha * (w)
    z = vcycle(z, r, A_3d, n);
    B = (r' * z) / (r_old' * z_old);
    p = z + B * p;   % Updating search
    w = A_3d*p;

    % er = ue - u
    % err_pcg(k) = norm(er)/norm(ue);

end;

% semilogy(erel_pcg,'b-',ms,9);
% legend('Jacobi', 'PCG');
