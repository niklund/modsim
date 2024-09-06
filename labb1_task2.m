
a = 1; b = 2;
fun = @(t, x) [1 - (b + 1)*x(1) + a*x(1)^(2)*x(2); b*x(1) - a*x(1)^(2)*x(2)];
x0 = 1; y0 = 1;
tspan = [0,1000];

rel_tol = 1e-6;
abs_tol = 1e-4;
opts = odeset('RelTol',rel_tol, 'AbsTol',abs_tol);
[t_out, y_out] = ode45(fun, tspan, [x0, y0], opts);

figure(1)
plot(t_out, y_out);
xlabel('Time (t)', interpreter = 'latex')
ylabel('Quantity of chemical (n)', interpreter = 'latex')
legend('Chemical A', 'Chemical B', interpreter = 'latex')

[x,y] = meshgrid(0:.08:3, 0:.08:3);

xdot = 1 - (b + 1).*x + a.*(x.^2).*y;
ydot = b.*x - a.*(x.^2).*y;

xdot = xdot./sqrt(xdot.^2 + ydot.^2);
ydot = ydot./sqrt(xdot.^2 + ydot.^2);

figure(2)
hold on
quiver(x,y, xdot, ydot)

contour(x, y, xdot, [0 0], 'r-');
contour(x, y, ydot, [0 0], 'g-');
legend('Quiver plot', 'xcline', 'ycline')

[peaks, index] = findpeaks(y_out(:,1));

T = (t_out(index(end)) - t_out(index(end-1)));










