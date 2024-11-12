hdr;

n=3; h=1/(n+1); h2i=1./(h*h);

x = h*[1:n]';
e = ones(n^3,1);
A = spdiags([-e 2*e -e], -1:1, n, n);

A_2d = kron(A, eye(n)) + kron(eye(n), A);
A_3d = kron(A, eye(n^2)) + kron(eye(n), A_2d);
% disp(full(A_3d))

A_3d = h2i * A_3d;

di = 6 / h^2;

k = [1:n]';
V = sqrt(2*h)*sin((h*pi)*(k*k'));
V_3d = kron(V, kron(V, V));

Lam = (2*h2i)*(1-cos(h*pi*k));
lmax = 2;
lmin = 0.6;

di = h*h/2;
LS = di*Lam;


nc = n/2;
xc = [1:nc]'/(nc+1);
J = lin_interp_mat([0; x; 1],[0; xc; 1]);
J = J(2:end-1,:);
% Aci = inv(J'*A*J);

J_3d = kron(J, kron(J, J));
Aci_3d = inv(J_3d'*A_3d*J_3d);

% ue=0*x;
ue=rand(n^3,1);
ue=V_3d*e;


b = A_3d*ue;
u=0*b; nsmooth=4; omega=2/3;

r=b;
for k=1:20
   e=ue-u;

%  u=u+cheb_smooth(r,lmax,lmin,nsmooth,A);
   for nu=1:nsmooth;
      u=u+omega*di*r;
      r=b-A_3d*u;
   end;

   % coarse grid correction
   % ec = J_3d * (Aci_3d * (J_3d' * r));
   % u = u + ec;
   % r = b - A_3d*u;

   err = norm(e)/norm(ue)
   if err <= 1e-8
     k
   end;


end;