n = 20;
n_iter = 20;
ue = rand(n^3,1);


y1 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 1);
y2 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 2);
y3 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 4);
y4 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 8);
x = 1:n_iter;

semilogy(x, y1, x, y2, x, y3, x, y4)
legend({'aspect ratio = 1mg','2mg', '3mg', '4mg'},'Location','southwest')




% y11 = ksp_mg_aniso_plot(n, n_iter, ue, 1, 1, 1);
% y21 = ksp_mg_aniso_plot(n, n_iter, ue, 1, 1, 2);
% y31 = ksp_mg_aniso_plot(n, n_iter, ue, 1, 1, 4);
%y41 = ksp_mg_aniso_plot(n, n_iter, ue, 1, 1, 8);
% x1 = 1:n_iter;

% semilogy(x1, y11, x1, y21, x1, y31, x1, y41)

% legend({'1ksp','2ksp', '4ksp', '8ksp'},'Location','southwest')
