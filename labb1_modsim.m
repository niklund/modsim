clear all

a=0.4; b=0.001; c=0.001; d=0.9; e=0.001; j=0.001; k=0.001; l=0.001; m=0.9; 
R0_1 = 600; R0_2 = 600; R0_3 = 100;
F0_1 = 400; F0_2 = 700; F0_3 = 100;
B0 = 10;
alpha=d/a; r0 = R0_1/(d/c); f0 = F0_1/(a/b);
r_grass = 5 * f0; r_sat = 5 * f0;

y0 = [R0_1, R0_2, R0_3; F0_1, F0_2, F0_3];
y0_mod = [r0, f0];
t = [0,200];
tau = t/(1/a);

rel_tol = 1e-6;
abs_tol = 1e-4;
opts = odeset('RelTol',rel_tol, 'AbsTol',abs_tol);

f = @(t, x) [a*x(1) - b*(x(1)*x(2)); c*x(1)*x(2) - d*x(2)];
g = @(t, x) [x(1)*(1-x(1)/r_grass) - x(1)*x(2)/(1 + x(1)/r_sat); alpha*(x(1)*x(2)/(1 + x(1)/r_sat) - x(2))];

h = @(t, x) [a*x(1) - b*(x(1)*x(2)) - e*(x(1)*x(3)); 
             c*x(1)*x(2) - d*x(2) - j*x(2)*x(3); 
             k*x(1)*x(3) + l*x(2)*x(3) - m*x(3)];

[T_out1, Y_out1] = ode45(f,t,y0(:,1),opts);
[T_out2, Y_out2] = ode45(f,t,y0(:,2),opts);
[T_out3, Y_out3] = ode45(f,t,y0(:,3),opts);

[T_out4, Y_out4] = ode45(g, tau, y0_mod, opts);
[T_out5, Y_out5] = ode45(h, t, [R0_1, F0_1, B0]', opts);

figure(1)
hold on
plot(T_out1, Y_out1);

xlabel('Time (t)', Interpreter='latex')
ylabel('Populations', Interpreter='latex')
legend('Rabbit pop.', 'Fox pop.', Interpreter='latex')
hold off
figure(2)
hold on
xlabel('Rabbit pop.', Interpreter='latex')
ylabel('Fox pop.', Interpreter='latex')
plot(Y_out1(:,1), Y_out1(:,2));
% plot(Y_out2(:,1), Y_out2(:,2));
% plot(Y_out3(:,1), Y_out3(:,2));
% legend('$R_0=600,F_0=400$', '$R_0=600,F_0=700$', '$R_0=100,F_0=100$', Interpreter='latex');
legend('Stochastic model', 'Continuous model', Interpreter='latex')
hold off

figure(3)
hold on
plot(T_out4, Y_out4);
xlabel('Time (t)', Interpreter='latex')
ylabel('Populations', Interpreter='latex')
legend('Rabbit pop.', 'Fox pop.', Interpreter='latex')
hold off

figure(4)
hold on
plot(Y_out4(:,1), Y_out4(:,2));
xlabel('Rabbit pop.', Interpreter='latex')
ylabel('Fox pop.', Interpreter='latex')
hold off

figure(5)
plot3(Y_out5(:,1), Y_out5(:,2), Y_out5(:,3));
xlabel('Rabbit pop.', Interpreter='latex')
ylabel('Fox pop.', Interpreter='latex')
zlabel('Bear pop.', Interpreter='latex')
