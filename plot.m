y1 = iso_multigrid_plot(4, 30);
y2 = aniso_multigrid_plot(4, 30, 1, 1, 2);
x = [1:30]

plot(x, y1, x, y2)