n = 20;
n_iter = 10;
ue = rand(n^3,1);

y1 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 1);
y2 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 2);
y3 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 4);
y4 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 8);
x = 1:n_iter;

semilogy(x, y1, x, y2, x, y3, x, y4)
legend({'1','2', '3', '4'},'Location','southwest')
