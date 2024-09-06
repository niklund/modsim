sigma = 10; beta = 8/3; rho = 50;
x0 = [-8, 8, 27]; tspan = [0 50];
delta = [10^(-15), 10^(-15), 10^(-15)];
pert = 10^(-15);

fun = @(t, x) [
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
expModel = fittype('pert * exp(lambda * x)', 'independent', 'x', 'coefficients', {'lambda'});
fitResult = fit(t_out, abs_diff, expModel, 'problem', pert);

est_lambda = fitResult.lambda;


figure(2)
hold on
plot(t_out, abs_diff)
plot(t_out, pert*exp(est_lambda*t_out))
hold off
