alpha = 77.27; beta = 8.375*10^(-6); gamma = 1.161;
tspan = [0, 1000];
x0 = [1, 2, 3];

fun = @(t, x) [
    alpha*(x(2) + x(1)*(1 - beta*x(1) - x(2)));
    (1/alpha)*(x(3) - (1 + x(1))*x(2));
    gamma*(x(1) - x(3))
];

syms x1 x2 x3;

sym_fun = [
    alpha*(x2 + x1*(1 - beta*x1 - x2));
    (1/alpha)*(x3 - (1 + x1)*x2);
    gamma*(x1 - x3)
];


Jac = jacobian(sym_fun, [x1, x2, x3]); 
Jac_fun = matlabFunction(Jac, 'Vars', {x1, x2, x3});
Jac = @(t, x) Jac_fun(x(1), x(2), x(3));


rel_tol = 1e-6;
abs_tol = 1e-4;
opts1 = odeset('RelTol',rel_tol, 'AbsTol',abs_tol);
opts2 = odeset('RelTol',rel_tol, 'AbsTol',abs_tol, Jacobian=Jac);
% 
% [t_out45, y_out45] = ode45(fun, tspan, x0, opts1);
% [t_out113, y_out113] = ode113(fun, tspan, x0, opts1);
[t_out23s, y_out23s] = ode23s(fun, tspan, x0, opts2);
[t_out15s, y_out15s] = ode23s(fun, tspan, x0, opts2);

figure(1)
plot3(y_out23s(:,1), y_out23s(:,2), y_out23s(:,3), LineWidth=2);
xlabel('Concentration $x_1$', Interpreter='latex')
ylabel('Concentration $x_2$', Interpreter='latex')
zlabel('Concentration $x_3$', Interpreter='latex')

figure(2)
hold on
yyaxis left;

plot(t_out15s, y_out15s(:,1))
plot(t_out15s, y_out15s(:,3))
yyaxis right;
plot(t_out15s, y_out15s(:,2))
legend('$x_1$', '$x_3$', '$x_2$', Interpreter='latex')
hold off

