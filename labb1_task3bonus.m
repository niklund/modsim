m1 = 1; m2 = 1; mu = m2/m1;
m = 1;
l = 1;
g = 9.82;
tspan = [0, 10];
X0 = [pi 0 1 1];
delta = [10^(-6) 0 0 0];

ode_dXdt = @(t, X) diffeq(X, m, l, g);

rel_tol = 1e-6;
abs_tol = 1e-4;
opts = odeset('RelTol',rel_tol, 'AbsTol',abs_tol);
[t_out, y_out] = ode45(ode_dXdt, tspan, X0, opts);
[t_out1, y_out1] = ode45(ode_dXdt, tspan, X0 - delta, opts);
[t_out2, y_out2] = ode45(ode_dXdt, tspan, X0 + delta, opts);

[x01, x02, y01, y02] = to_cartesian(y_out, l);
[x11, x12, y11, y12] = to_cartesian(y_out1, l);
[x21, x22, y21, y22] = to_cartesian(y_out2, l);


figure(1)
hold on
plot(x02, y02)
plot(x12, y12)
plot(x22, y22)
xlabel('x', Interpreter='latex')
ylabel('y', Interpreter='latex')
legend('with $\theta_0$', 'with $\theta_0 - \delta$', 'with $\theta_0 + \delta$', Interpreter='latex')
hold off


function dXdt = diffeq(X, m, l, g)
    
    denom = m*l^2*(16 - 9*cos(X(1) - X(2))^2);
    prefac = 6/(m*l^2);
    thetadot1 = (prefac/denom)*(2*X(3) - 3*cos(X(1) - X(2))*X(4));
    thetadot2 = (prefac/denom)*(8*X(4) - 3*cos(X(1) - X(2))*X(3));

    dXdt = [
        thetadot1;
        thetadot2;
        -0.5*m*l^2*(thetadot1*thetadot2*sin(X(1) - X(2)) + 3*(g/l)*sin(X(1)));
        -0.5*m*l^2*(-thetadot1*thetadot2*sin(X(1) - X(2)) + (g/l)*sin(X(2)));
    ];

end

function [x1, x2, y1, y2] = to_cartesian(sol, l)
    x1 = l*sin(sol(:,1));
    x2 = x1 + l*sin(sol(:,2));
    y1 = -l*cos(sol(:,1));
    y2 = y1 - l*cos(sol(:,2));
end
    