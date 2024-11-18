n = 15;
n_iter = 50;
ue = rand(n^3,1);

y1 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 1);
y2 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 2);
y3 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 4);
y4 = aniso_multigrid_plot(n, n_iter, ue, 1, 1, 8);
x = 1:n_iter;


semilogy(x, y1, 'Linewidth', 2,  x, y2, 'Linewidth', 2, x, y3, 'Linewidth', 2, x, y4, 'Linewidth', 2);
legend({'r=1','r=2', 'r=4', 'r=8'},'Location','southwest')
ylabel('norm of residual');
xlabel('number of iterations')
