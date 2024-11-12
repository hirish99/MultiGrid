hdr;

n=200; h=1/(n+1); h2i=1./(h*h);

x = h*[1:n]';
e = ones(n,1);
A = h2i*spdiags([-e 2*e -e], -1:1, n, n);

k = [1:n]';
V = sqrt(2*h)*sin((h*pi)*(k*k'));
Lam = (2*h2i)*(1-cos(h*pi*k));
lmax = 2;
lmin = 0.6;

di = h*h/2;
LS = di*Lam;

ue=0*x;
ue=rand(n,1);
ue=V*e;


b = A*ue; u=0*b; nsmooth=4; omega=2/3;

r=b;
for k=1:20
   e=ue-u; eh=V'*e; uh=V'*u;
   %plot(LS,0*LS,'k-',LS,eh,'r-',lw,3);
   %title(['Jac: $k$ = ' num2str([k-1 n]) ', $\omega$ = ' num2str(omega)],intp,ltx,fs,22);
   %xlabel(['$\lambda(D^{-1}A)$'],intp,ltx,fs,18);
   %ylabel(['${\hat e_k}/\|{\underline u}\|$'],intp,ltx,fs,18);

   %pause(.05 + .95/k);
   %plot(LS,0*LS,'k-',LS,eh,'b-',lw,3); hold on;

%  u=u+cheb_smooth(r,lmax,lmin,nsmooth,A);
   for nu=1:nsmooth;
      u=u+omega*di*r;
      r=b-A*u;
   end;

   % coarse grid correction
   xc = 2/(201) * [1:100]';
   J = lin_interp_mat(x, xc);
   Ac = J' * A * J;
   ec = J * inv(Ac) * J' * r;
   u = u + ec;
   r = b - A*u;

   if norm(e)/norm(ue) <= 1e-8
     k
   endif

end;

plot(LS,0*LS,'k-',LS,eh,'r-',lw,3); hold on;
title(['Jac: $k$ = ' num2str([k-1 n]) ', $\omega$ = ' num2str(omega)],intp,ltx,fs,22);
xlabel(['$\lambda(D^{-1}A)$'],intp,ltx,fs,18);
ylabel(['${\hat e_k}/\|{\underline u}\|$'],intp,ltx,fs,18);
