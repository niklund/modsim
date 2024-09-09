sigma = 10; beta = 8/3; rho = 50;
x0 = [-8, 8, 27]; tspan = [0 50];
delta = [10^(-15), 10^(-15), 10^(-15)];
pert = 10^(-15);

fun = @(t, x) [
    sigma*(x(2) - x(1));
    x(1)*(rho - x(3)) - x(2);
    x(1)*x(2) - beta*x(3);
];

fun2 = @(x) [
    sigma*(x(2) - x(1));
    x(1)*(rho - x(3)) - x(2);
    x(1)*x(2) - beta*x(3);
];

rel_tol = 1e-6;
abs_tol = 1e-4;
opts = odeset('RelTol',rel_tol, 'AbsTol',abs_tol);
[t_out, y_out] = ode45(fun, tspan, x0, opts);
[t_out2, y_out2] = ode45(fun, tspan, x0 + delta, opts);

figure(1)
hold on
plot3(y_out(:,1), y_out(:,2), y_out(:,3), color="#4DBEEE")
plot3(y_out2(:,1), y_out2(:,2), y_out2(:,3), color="#EDB120")
legend('$x_0$','$x_0 + \delta$', Interpreter='latex')
xlabel('$x_1$', Interpreter='latex')
ylabel('$x_2$', Interpreter='latex')
zlabel('$x_3$', Interpreter='latex')
hold off
view(3)

len_y_out = length(y_out(:,1));
abs_diff = vecnorm(y_out2(1:len_y_out,:) - y_out, 2, 2);

logdiff = log(abs_diff/pert);

mf = fit(t_out(1650:5100), logdiff(1650:5100), 'poly1');

figure(2)
hold on
plot(t_out, pert*exp(mf.p1*t_out))
hold off

X01 = [9, 15, 47]; X02 = [-9, -15, 47]; X03 = [6, 14, 22];
fopts = optimset('Display','off');
X1 = fsolve(fun2, X01, fopts);
X2 = fsolve(fun2, X02, fopts);

syms x1 x2 x3;

sym_fun = [
    sigma*(x2 - x1);
    x1*(rho - x3) - x2;
    x1*x2 - beta*x3;
];

Jac = jacobian(sym_fun, [x1, x2, x3]); 
Jac_fun = matlabFunction(Jac, 'Vars', {x1, x2, x3});
Jac = @(t, x) Jac_fun(x(1), x(2), x(3));

jac_eigs = eig(Jac(0, X1));


