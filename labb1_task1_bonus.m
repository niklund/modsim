
a=0.4; b=0.001; c=0.001; d=0.9;
R = 600;
F = 400;
t = 0;

iters = 10^5;
R_vec = zeros(1,iters);
F_vec = zeros(1,iters);
R_vec(1) = R;
F_vec(1) = F;
t_vec = zeros(1, iters);

for i = 2:iters
    k1 = a*R;
    k2 = b*F*R;
    k3 = d*F;
    k_tot = k1 + k2 + k3;
    
    r0 = rand();
    dt = -(1/k_tot)*log(r0);
    r1 = rand();

    if r1 < k1/k_tot
        R = R + 1;
    
        elseif (r1 < (k1 + k2)/k_tot)
            R = R - 1;
            F = F + 1;
    
        else 
            F = F - 1;
    end
    R_vec(i) = R;
    F_vec(i) = F;
    t_vec(i) = t_vec(i-1) + dt;
end

figure(1)
hold on
xlabel('Time (t)', Interpreter='latex')
ylabel('Populations', Interpreter='latex')
plot(t_vec, R_vec, t_vec, F_vec);
legend('Rabbit pop.', 'Fox pop.', Interpreter='latex')
hold off
figure(2)
hold on
plot(R_vec, F_vec);
xlabel('Rabbit pop.', Interpreter='latex')
ylabel('Fox pop.', Interpreter='latex')
hold off

biters = 10^5;
Rc_vec = zeros(10, biters);
Rc_vec(:,1) = 600;
Fc_vec = zeros(10, biters);
Fc_vec(:,1) = [2;4;6;8;10;12;14;16;18;20];

probvec = zeros(50,10);
for k = 1:50
for i = 1:10
    Rc = Rc_vec(i,1);
    Fc = Fc_vec(i,1);
    for j = 2:biters
        k1 = a*Rc;
        k2 = b*Fc*Rc;
        k3 = d*Fc;
        k_tot = k1 + k2 + k3;
    
        r0 = rand();
        dt = -(1/k_tot)*log(r0);
        r1 = rand();

        if r1 < k1/k_tot
            Rc = Rc + 1;
           
        elseif (r1 < (k1 + k2)/k_tot)
              Rc = Rc - 1;
              Fc = Fc + 1;
        
        else 
              Fc = Fc - 1;
        end
    Rc_vec(i, j) = Rc;
    Fc_vec(i, j) = Fc;
    end
    if Fc_vec(i, end) == 0
        probvec(k, i) = 1;
    end
end
end

prob = sum(probvec)/50;
figure(3)
plot(Fc_vec(:,1), prob)
xlabel('Initial fox pop.', Interpreter='latex')
ylabel('Probability of fox pop dying', Interpreter='latex')
