%% SC42056 Optimization for Systems and Control
% Nonlinear Programming Assignment - Group 22
% Sven Rutgers - 4600150
% Melis Orhun  - 4912071
% =============================================================
% Optimization algorithm with continuous and discrete r
% *** RudiMINtal ****
% =============================================================
% *** TABLE OF CONTENTS ***
%  
% 1) PARAMETERS
% 2) OPTIMIZATION - Continuous r(k)
%  * fmincon with SQP and r=0
%  * fmincon with Interior-point and r=0
%  * fmincon with SQP and r=0.99
%  * fmincon with Interior-point and r=0.99
%  * fmincon with SQP scaled
% 3) OPTIMIZATION - Discrete r(k)
%  * genetic algorithm
%  * no control case calculation
% 3) PLOTS
%  * TTS, Elapsed Time and Num. of Iterations Table
%  * Plots
% 5) FUNCTIONS FOR CONTINUOUS r(k)
%  * nonlinear constraints function definition
%  * cost function definition
% 6) FUNCTIONS FOR DISCRETE r(k)
%  * nonlinear constraints function definition
%  * cost function definition
% 7) FUNCTIONS FOR PLOTTING
%  * output optimum states function definition
% =============================================================

close all
clear
clc

% =============================================================
% 1) PARAMETERS
% =============================================================

% Group-specific parameters
E1 = 0 + 1;
E2 = 7 + 5;

% Model Parameters
Li = 1;              % km
T = 10/3600;         % h
tau = 18/3600;       % h
mu = 60;             % km^2/h
Cr = 2000;           % veh/h
K = 0.04;            % veh/(km*lane)
v_f = 100;           % km/h
rho_m = 180;         % veh/(km*lane)
Dr = 1500;           % veh/h
a = 1.87;
lmbda = 4;
rho_c = 33.5 + E1/3;

% =============================================================
% 2) OPTIMIZATION - Continuous r(k)
% =============================================================

% Define upper and lower bounds
lb_fmincon = zeros(121,1);
ub_fmincon = ones(121,1);

% Set initial points
r0 = 0 * ones(121,1);

% Run optimization algorithms
% SQP with r = 0
options_fmincon = optimoptions('fmincon','Display','iter','Algorithm','sqp','PlotFcn','optimplotfval'); % Use sqp
t_begin_sqp_0 = tic;
[r_sqp_0, out_sqp_0, exitflag_sqp_0, data_sqp_0] = fmincon(@cost_function_fmincon,r0,[],[],[],[],lb_fmincon,ub_fmincon,@nlcon_fmincon,options_fmincon);
t_sqp_0 = toc(t_begin_sqp_0);
[rho_sqp_0, v_sqp_0, wr_sqp_0] = optimum_states(r_sqp_0);

% Interior Point with r = 0
options_fmincon = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',inf,'MaxIterations',inf,'PlotFcn','optimplotfval'); 
t_begin_ip_0 = tic;
[r_ip_0, out_ip_0, exitflag_ip_0, data_ip_0] = fmincon(@cost_function_fmincon,r0,[],[],[],[],lb_fmincon,ub_fmincon,@nlcon_fmincon,options_fmincon);
t_ip_0 = toc(t_begin_ip_0);
[rho_ip_0, v_ip_0, wr_ip_0] = optimum_states(r_ip_0);

% Set initial point r = 0.99
r0 = 0.99 * ones(121,1);

% SQP with r = 0.99
options_fmincon = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'StepTolerance',1e-90000000000,'PlotFcn','optimplotfval'); % Use sqp
t_begin_sqp_099 = tic;
[r_sqp_099, out_sqp_099, exitflag_sqp_099, data_sqp_099] = fmincon(@cost_function_fmincon,r0,[],[],[],[],lb_fmincon,ub_fmincon,@nlcon_fmincon,options_fmincon);
t_sqp_099 = toc(t_begin_sqp_099);
[rho_sqp_099, v_sqp_099, wr_sqp_099] = optimum_states(r_sqp_099);

% Interior Point with r = 0
options_fmincon = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',inf,'MaxIterations',inf,'PlotFcn','optimplotfval'); 
t_begin_ip_099 = tic;
[r_ip_099, out_ip_099, exitflag_ip_099, data_ip_099] = fmincon(@cost_function_fmincon,r0,[],[],[],[],lb_fmincon,ub_fmincon,@nlcon_fmincon,options_fmincon);
t_ip_099 = toc(t_begin_ip_099);
[rho_ip_099, v_ip_099, wr_ip_099] = optimum_states(r_ip_099);

% SQP with r = 0 - scaled upper and lower bounds
lb_fmincon = 0.15*ones(121,1);
ub_fmincon = 0.45*ones(121,1);
options_fmincon = optimoptions('fmincon','Display','iter','Algorithm','sqp');
t_begin_sqp_0_scaled = tic;
[r_sqp_0_scaled, out_sqp_0_scaled, exitflag_sqp_0_scaled, data_sqp_0_scaled] = fmincon(@cost_function_fmincon,r0,[],[],[],[],lb_fmincon,ub_fmincon,@nlcon_fmincon,options_fmincon);
t_sqp_0_scaled = toc(t_begin_sqp_0_scaled);
[rho_sqp_0_scaled, v_sqp_0_scaled, wr_sqp_0_scaled] = optimum_states(r_sqp_0_scaled);

% =============================================================
% 3) OPTIMIZATION - Discrete r(k) 
% *** NOTE *** 
% If the exit flag for the genetic algorithm is negative
% re-run the code. Due to the heuristic nature of the algorithm, it can
% not always find an optimum solution.
% =============================================================
%
% Define upper and lower bounds
lb_ga = ones(121,1);
ub_ga = 7*ones(121,1);

% Run optimization algorithm
options_ga = optimoptions('ga','PlotFcn',@gaplotbestf,'Display','iter','ConstraintTolerance',1e-10, ...
                         'FunctionTolerance',1e-10,'MaxGenerations',1000);
t_begin_ga = tic;
[r_optimum_ga, out_ga, exitflag_ga, data_ga]=ga(@cost_function_ga,121,[],[],[],[],lb_ga,ub_ga,@nlcon_ga,1:1:121,options_ga);
t_ga = toc(t_begin_ga);

% Initialize the optimum r matrix
r_ga = zeros(1,121);
for i = 1:1:121
    r_ga(i) = mapVariable(r_optimum_ga(i));
end

[rho_ga, v_ga, wr_ga] = optimum_states(r_ga);


% No control case calculations
r = 1;
% Initialize the output matrices
rho_no_control = zeros(4,121);
v_no_control = zeros(4,121);
wr_no_control = zeros(1,121);
y = 0;

% Set initial conditions
rho_no_control(1:4,1) = [30; 30; 30; 30];
v_no_control(1:4,1) = [80; 80; 80; 80];
wr_no_control(1) = 0;
q0 = 8000 + 100*E1;

for j = 1:1:120
    
    state_update_input = [rho_no_control(1:4,j); v_no_control(1:4,j); wr_no_control(j); rho_c; r; q0];

    [rho_no_control(1:4,j+1), v_no_control(1:4,j+1), wr_no_control(j+1)] = state_update_no_control(state_update_input);

    if j == 30
        q0 = 4000 + 100*E2;
    end

    y = y + T*wr_no_control(j) + T*Li*lmbda*( rho_no_control(1,j) + rho_no_control(2,j) + rho_no_control(3,j) + rho_no_control(4,j) );
end
y = y + T*wr_no_control(j+1) + T*Li*lmbda*( rho_no_control(1,j+1) + rho_no_control(2,j+1) + rho_no_control(3,j+1) + rho_no_control(4,j+1) );

% =============================================================
% 4) PLOTS
% =============================================================

% Display data of different algorithms
Tbl = table([out_sqp_0;out_ip_0;out_sqp_099;out_ip_099;out_ga;y],[t_sqp_0;t_ip_0;t_sqp_099;t_ip_099;t_ga;0], ...
          [data_sqp_0.iterations;data_ip_0.iterations;data_sqp_099.iterations;data_ip_099.iterations;data_ga.generations;0], ...
          [exitflag_sqp_0;exitflag_ip_0;exitflag_sqp_099;exitflag_ip_099;exitflag_ga;0], ...
          'VariableNames',{'Mini TTS (h)','Elapsed Time (s)','Iterations/Generations','ExitFlag'}, ...
          'RowNames',{'RudiMINtal SQP (r=0)','RudiMINtal IP (r=0)','RudiMINtal SQP (r=0.99)','RudiMINtal IP (r=0.99)','RudiMINtal Genetic Algorithm','No Control Case'});
disp(Tbl)

set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex'); 
set(groot,'defaultLineLineWidth',1)
set(groot,'defaultAxesFontSize',11)
set(groot,'defaultAxesFontWeight',"normal")

% Plot the optimized solutions for RudiMINtal
figure(1) % SQP - density r=0
plot(0:10:1200,rho_sqp_0(1,1:121),0:10:1200,rho_sqp_0(2,1:121),0:10:1200,rho_sqp_0(3,1:121),0:10:1200,rho_sqp_0(4,1:121))
title("Traffic Density vs. Time for Each Lane (with SQP and $r_{sp}(k)=0$)")
xlabel("Time (s)"); ylabel("Traffic Density [veh/(km.lane)]"); grid minor
legend('$\rho_1$(t)','$\rho_2$(t)','$\rho_3$(t)','$\rho_4$(t)')

figure(2) % SQP - mean speed r=0
plot(0:10:1200,v_sqp_0(1,1:121),0:10:1200,v_sqp_0(2,1:121),0:10:1200,v_sqp_0(3,1:121),0:10:1200,v_sqp_0(4,1:121))
title("Mean Speed vs. Time (with SQP and $r_{sp}(k)=0$)")
xlabel("Time (s)"); ylabel("Mean Speed [km/h]"); grid minor
legend('$v_1$(t)','$v_2$(t)','$v_3$(t)','$v_4$(t)','Location','northwest')

figure(3) % SQP - number of vehicles in queue on the on-ramp r=0
plot(0:10:1200,wr_sqp_0) 
title(sprintf("Number of Vehicles in Queue on the On-ramp vs. Time\n(with SQP and $r_{sp}(k)=0$)"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor

figure(4) % SQP - ramp metering rate r=0
plot(0:10:1200,r_sqp_0)
title("Ramp Metering Rate vs. Time (with SQP and $r_{sp}(k)=0$)")
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor

figure(5) % IP - density r=0.99
plot(0:10:1200,rho_ip_099(1,1:121),0:10:1200,rho_ip_099(2,1:121),0:10:1200,rho_ip_099(3,1:121),0:10:1200,rho_ip_099(4,1:121))
title("Traffic Density vs. Time for Each Lane (IP with $r_{sp}(k)=0.99$)")
xlabel("Time (s)"); ylabel("Traffic Density [veh/(km.lane)]"); grid minor
legend('$\rho_1$(t)','$\rho_2$(t)','$\rho_3$(t)','$\rho_4$(t)')

figure(6) % IP - mean speed r=0.99
plot(0:10:1200,v_ip_099(1,1:121),0:10:1200,v_ip_099(2,1:121),0:10:1200,v_ip_099(3,1:121),0:10:1200,v_ip_099(4,1:121))
title("Mean Speed vs. Time (IP with $r_{sp}(k)=0.99$)")
xlabel("Time (s)"); ylabel("Mean Speed [km/h]"); grid minor
legend('$v_1$(t)','$v_2$(t)','$v_3$(t)','$v_4$(t)','Location','northwest')

figure(7) % IP - number of vehicles in queue on the on-ramp r=0.99
plot(0:10:1200,wr_ip_099) 
title(sprintf("Number of Vehicles in Queue on the On-ramp vs. Time\n(IP with $r_{sp}(k)=0.99$)"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor

figure(8) % IP - ramp metering rate r=0.99
plot(0:10:1200,r_ip_099)
title("Ramp Metering Rate vs. Time (IP with $r_{sp}(k)=0.99$)")
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor

% Plot the optimized solutions for SQP - IP comparison r=0
figure(9) % density comparison
sgtitle(sprintf("Traffic Density Comparison\nfor SQP and IP Solutions of RudiMINtal with $r_{sp}(k)=0$"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,rho_ip_0(1,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(1,1:121),":");
title("Traffic Density vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{1,ip}$(t)','$\rho_{1,sqp}$(t)')
%
subplot(2,2,2)
plot(0:10:1200,rho_ip_0(2,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(2,1:121),":");
title("Traffic Density vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{2,ip}$(t)','$\rho_{2,sqp}$(t)')
%
subplot(2,2,3)
plot(0:10:1200,rho_ip_0(3,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(3,1:121),":");
title("Traffic Density vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{3,ip}$(t)','$\rho_{3,sqp}$(t)')
%
subplot(2,2,4)
plot(0:10:1200,rho_ip_0(4,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(4,1:121),":");
title("Traffic Density vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{4,ip}$(t)','$\rho_{4,sqp}$(t)','Location','southwest')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(10) % speed comparison
sgtitle(sprintf("Mean Speed Comparison\nfor SQP and IP Solutions of RudiMINtal with $r_{sp}(k)=0$"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,v_ip_0(1,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(1,1:121),":");
title("Mean Speed vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{1,ip}$(t)','$v_{1,sqp}$(t)','Location','southeast')
%
subplot(2,2,2)
plot(0:10:1200,v_ip_0(2,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(2,1:121),":");
title("Mean Speed vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{2,ip}$(t)','$v_{2,sqp}$(t)','Location','southeast')
%
subplot(2,2,3)
plot(0:10:1200,v_ip_0(3,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(3,1:121),":");
title("Mean Speed vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{3,ip}$(t)','$v_{3,sqp}$(t)')
%
subplot(2,2,4)
plot(0:10:1200,v_ip_0(4,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(4,1:121),":");
title("Mean Speed vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{4,ip}$(t)','$v_{4,sqp}$(t)')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(11) % queue length comparison
plot(0:10:1200,wr_ip_0,"--"); hold on; plot(0:10:1200,wr_sqp_0,":"); 
title(sprintf("Queue Length vs. Time\nfor SQP and IP Solutions of RudiMINtal with $r_{sp}(k)=0$"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor
legend('$w_{r,ip}$(t)','$w_{r,sqp}$(t)')

figure(12) % ramp metering rate comparison
plot(0:10:1200,r_ip_0,"--"); hold on; plot(0:10:1200,r_sqp_0,":"); 
title(sprintf("Ramp Metering Rate vs. Time\n for SQP and IP Solutions of RudiMINtal with $r_{sp}(k)=0$"))
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor
legend('$r_{ip}$(t)','$r_{sqp}$(t)','Location','southeast')

% Plot the optimized and no control solutions
figure(13) % density comparison
sgtitle(sprintf("Traffic Density Comparison\nfor No-Control and SQP-Optimized Solutions of RudiMINtal"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,rho_no_control(1,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(1,1:121),":");
title("Traffic Density vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
%
subplot(2,2,2)
plot(0:10:1200,rho_no_control(2,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(2,1:121),":");
title("Traffic Density vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
%
subplot(2,2,3)
plot(0:10:1200,rho_no_control(3,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(3,1:121),":");
title("Traffic Density vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
%
subplot(2,2,4)
plot(0:10:1200,rho_no_control(4,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(4,1:121),":");
title("Traffic Density vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(14) % speed comparison
sgtitle(sprintf("Mean Speed Comparison\nfor No-Control and SQP-Optimized Solutions of RudiMINtal"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,v_no_control(1,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(1,1:121),":");
title("Mean Speed vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','southwest')
%
subplot(2,2,2)
plot(0:10:1200,v_no_control(2,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(2,1:121),":");
title("Mean Speed vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','southwest')
%
subplot(2,2,3)
plot(0:10:1200,v_no_control(3,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(3,1:121),":");
title("Mean Speed vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','southwest')
%
subplot(2,2,4)
plot(0:10:1200,v_no_control(4,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(4,1:121),":");
title("Mean Speed vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','east')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(15) % queue length comparison
plot(0:10:1200,wr_no_control,"--"); hold on; plot(0:10:1200,wr_sqp_0,":"); 
title(sprintf("Queue Length vs. Time\nfor No-Control and SQP-Optimized Solutions of RudiMINtal"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$')

% Plot the SQP, NoControl and genetic algorithm solutions
figure(16) % density comparison
sgtitle(sprintf("Traffic Density Comparison for No-Control,\nSQP-Optimized and GA-Optimized Solutions of RudiMINtal"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,rho_no_control(1,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(1,1:121),":"); plot(0:10:1200,rho_ga(1,1:121));
title("Traffic Density vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','northwest')
%
subplot(2,2,2)
plot(0:10:1200,rho_no_control(2,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(2,1:121),":"); plot(0:10:1200,rho_ga(2,1:121));
title("Traffic Density vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','northwest')
%
subplot(2,2,3)
plot(0:10:1200,rho_no_control(3,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(3,1:121),":"); plot(0:10:1200,rho_ga(3,1:121));
title("Traffic Density vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','northwest')
%
subplot(2,2,4)
plot(0:10:1200,rho_no_control(4,1:121),"--"); hold on; plot(0:10:1200,rho_sqp_0(4,1:121),":"); plot(0:10:1200,rho_ga(4,1:121));
title("Traffic Density vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','northwest')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(17) % speed comparison
sgtitle(sprintf("Mean Speed Comparison for No-Control,\nSQP-Optimized and GA-Optimized Solutions of RudiMINtal"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,v_no_control(1,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(1,1:121),":"); plot(0:10:1200,v_ga(1,1:121));
title("Mean Speed vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','southwest')
%
subplot(2,2,2)
plot(0:10:1200,v_no_control(2,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(2,1:121),":"); plot(0:10:1200,v_ga(1,1:121));
title("Mean Speed vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','southwest')
%
subplot(2,2,3)
plot(0:10:1200,v_no_control(3,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(3,1:121),":"); plot(0:10:1200,v_ga(1,1:121));
title("Mean Speed vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','southwest')
%
subplot(2,2,4)
plot(0:10:1200,v_no_control(4,1:121),"--"); hold on; plot(0:10:1200,v_sqp_0(4,1:121),":"); plot(0:10:1200,v_ga(1,1:121));
title("Mean Speed vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','southwest')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(18) % queue length comparison
plot(0:10:1200,wr_no_control,"--"); hold on; plot(0:10:1200,wr_sqp_0,":"); plot(0:10:1200,wr_ga);
title(sprintf("Queue Length vs. Time for No-Control,\nSQP-Optimized and GA-Optimized Solutions of RudiMINtal"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA')

figure(19) % ramp metering rate comparison
plot(0:10:1200,ones(1,121),"--"); hold on; plot(0:10:1200,r_sqp_0,":"); scatter(0:10:1200,r_ga); 
title(sprintf("Ramp Metering Rate vs. Time for No-Control,\nSQP-Optimized and GA-Optimized Solutions of RudiMINtal"))
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor
legend('No Control','SQP ($r_{sp}(k)=0$)','GA','Location','southeast')

figure(20) % ramp metering rate comparison
hold on; plot(0:10:1200,r_sqp_0_scaled,":"); scatter(0:10:1200,r_ga); 
title(sprintf("Ramp Metering Rate vs. Time for \nSQP-Optimized (Scaled) and GA-Optimized Solutions of RudiMINtal"))
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor
legend('SQP ($r_{sp}(k)=0$)','GA','Location','east')

% =============================================================
% 5) FUNCTIONS FOR CONTINUOUS r(k)
% =============================================================

% Nonlinear constraints function
function [c, ceq] = nlcon_fmincon(r)

    % Group-specific parameters
    E1 = 0 + 1;
    E2 = 7 + 5;

    % Model Parameters
    Li = 1;              % km
    T = 10/3600;         % h
    tau = 18/3600;       % h
    mu = 60;             % km^2/h
    Cr = 2000;           % veh/h
    K = 0.04;            % veh/(km*lane)
    v_f = 100;           % km/h
    rho_m = 180;         % veh/(km*lane)
    Dr = 1500;           % veh/h
    a = 1.87;
    lmbda = 4;
    rho_c = 33.5 + E1/3; % veh/(km*lane)

    % Initialize the output matrices
    rho = zeros(4,121);
    v = zeros(4,121);
    wr = zeros(1,121);
    q = zeros(4,1);
    V = zeros(4,1);
    c = zeros(484,1); % Nonlinear inequality constraints
    ceq = []; % Nonlinear equality constraints

    % Initial input
    q0 = 8000 + 100*E1;

    % Initial states
    rho(1:4,1) = 30*ones(4,1);
    v(1:4,1) = 80*ones(4,1);

    for j = 1:1:121
        
        % Change input at k=30
        if (j >= 31)
            q0 = 4000 + 100*E2;
        end
        
        % Calculate traffic flow for each segment
        q(1) = lmbda*(rho(1,j)*v(1,j));
        q(2) = lmbda*(rho(2,j)*v(2,j));
        q(3) = lmbda*(rho(3,j)*v(3,j));
        q(4) = lmbda*(rho(4,j)*v(4,j));
    
        % Calculate the desired speed for each segment
        V(1) = v_f * exp( (-1/a) * ( rho(1,j)/rho_c )^a );
        V(2) = v_f * exp( (-1/a) * ( rho(2,j)/rho_c )^a );
        V(3) = v_f * exp( (-1/a) * ( rho(3,j)/rho_c )^a );
        V(4) = v_f * exp( (-1/a) * ( rho(4,j)/rho_c )^a );

        % Calculate the traffic flow that enters segment 4 from the on-ramp
        qr4 = min( [r(j)*Cr, Dr+(wr(j)/T), Cr*(rho_m-rho(4,j))/(rho_m-rho_c)] );
    
        % Calculate the state update for density
        rho(1,j+1) = rho(1,j) + (T/(lmbda*Li)) * (q0-q(1));
        rho(2,j+1) = rho(2,j) + (T/(lmbda*Li)) * (q(1)-q(2));
        rho(3,j+1) = rho(3,j) + (T/(lmbda*Li)) * (q(2)-q(3));
        rho(4,j+1) = rho(4,j) + (T/(lmbda*Li)) * (q(3)-q(4)+qr4);
        
        % Calculate the state update for velocity
        v(1,j+1) = v(1,j) + (T/tau) * (V(1)-v(1,j)) - (mu*T/(tau*Li)) * (rho(2,j)-rho(1,j)) / (rho(1,j)+K);
        v(2,j+1) = v(2,j) + (T/tau) * (V(2)-v(2,j)) + (T/Li) * v(2,j) * (v(1,j)-v(2,j)) - (mu*T/(tau*Li)) * (rho(3,j)-rho(2,j)) / (rho(2,j)+K);
        v(3,j+1) = v(3,j) + (T/tau) * (V(3)-v(3,j)) + (T/Li) * v(3,j) * (v(2,j)-v(3,j)) - (mu*T/(tau*Li)) * (rho(4,j)-rho(3,j)) / (rho(3,j)+K);
        v(4,j+1) = v(4,j) + (T/tau) * (V(4)-v(4,j)) + (T/Li) * v(4,j) * (v(3,j)-v(4,j));

        % Calculate the state update for queue length
        wr(j+1) = wr(j) + T*(Dr-qr4);

        c((j-1)*4+1,1) = -rho_c + rho(1,j+1);
        c((j-1)*4+2,1) = -rho_c + rho(2,j+1);
        c((j-1)*4+3,1) = -rho_c + rho(3,j+1);
        c((j-1)*4+4,1) = -rho_c + rho(4,j+1);

    end
    
end

% Cost function
function y = cost_function_fmincon(r)

    % Group-specific parameters
    E1 = 0 + 1;
    E2 = 7 + 5;

    % Model Parameters
    Li = 1;              % km
    T = 10/3600;         % h
    tau = 18/3600;       % h
    mu = 60;             % km^2/h
    Cr = 2000;           % veh/h
    K = 0.04;            % veh/(km*lane)
    v_f = 100;           % km/h
    rho_m = 180;         % veh/(km*lane)
    Dr = 1500;           % veh/h
    a = 1.87;
    lmbda = 4;
    rho_c = 33.5 + E1/3; % veh/(km*lane)

    % Initialize the output matrices
    rho = zeros(4,122);
    v = zeros(4,122);
    wr = zeros(1,122);
    q = zeros(4,1);
    V = zeros(4,1);

    % Initial input
    q0 = 8000 + 100*E1;

    % Initial states
    rho(1:4,1) = 30*ones(4,1);
    v(1:4,1) = 80*ones(4,1);

    % Initialize output
    y = 0;

    for j = 1:1:121

        % Change input at k=30
        if (j >= 31)
            q0 = 4000 + 100*E2;
        end
        
        % Calculate traffic flow for each segment
        q(1) = lmbda*(rho(1,j)*v(1,j));
        q(2) = lmbda*(rho(2,j)*v(2,j));
        q(3) = lmbda*(rho(3,j)*v(3,j));
        q(4) = lmbda*(rho(4,j)*v(4,j));
    
        % Calculate the desired speed for each segment
        V(1) = v_f * exp( (-1/a) * ( rho(1,j)/rho_c )^a );
        V(2) = v_f * exp( (-1/a) * ( rho(2,j)/rho_c )^a );
        V(3) = v_f * exp( (-1/a) * ( rho(3,j)/rho_c )^a );
        V(4) = v_f * exp( (-1/a) * ( rho(4,j)/rho_c )^a );
        
        % Calculate the traffic flow that enters segment 4 from the on-ramp
        qr4 = min( [r(j)*Cr, Dr+(wr(j)/T), Cr*(rho_m-rho(4,j))/(rho_m-rho_c)] );
    
        % Calculate the state update for density
        rho(1,j+1) = rho(1,j) + (T/(lmbda*Li)) * (q0-q(1));
        rho(2,j+1) = rho(2,j) + (T/(lmbda*Li)) * (q(1)-q(2));
        rho(3,j+1) = rho(3,j) + (T/(lmbda*Li)) * (q(2)-q(3));
        rho(4,j+1) = rho(4,j) + (T/(lmbda*Li)) * (q(3)-q(4)+qr4);
        
        % Calculate the state update for velocity
        v(1,j+1) = v(1,j) + (T/tau) * (V(1)-v(1,j)) - (mu*T/(tau*Li)) * (rho(2,j)-rho(1,j)) / (rho(1,j)+K);
        v(2,j+1) = v(2,j) + (T/tau) * (V(2)-v(2,j)) + (T/Li) * v(2,j) * (v(1,j)-v(2,j)) - (mu*T/(tau*Li)) * (rho(3,j)-rho(2,j)) / (rho(2,j)+K);
        v(3,j+1) = v(3,j) + (T/tau) * (V(3)-v(3,j)) + (T/Li) * v(3,j) * (v(2,j)-v(3,j)) - (mu*T/(tau*Li)) * (rho(4,j)-rho(3,j)) / (rho(3,j)+K);
        v(4,j+1) = v(4,j) + (T/tau) * (V(4)-v(4,j)) + (T/Li) * v(4,j) * (v(3,j)-v(4,j));
    
        % Calculate the state update for queue length
        wr(j+1) = wr(j) + T * (Dr-qr4);
        
        y = y + T*Li*lmbda*(rho(1,j) + rho(2,j) + rho(3,j) + rho(4,j)) + T*wr(j);

    end

end

% =============================================================
% 6) FUNCTIONS FOR DISCRETE r(k)
% =============================================================

% Nonlinear constraints function
function [c, ceq] = nlcon_ga(r)

    % Group-specific parameters
    E1 = 0 + 1;
    E2 = 7 + 5;

    % Model Parameters
    Li = 1;              % km
    T = 10/3600;         % h
    tau = 18/3600;       % h
    mu = 60;             % km^2/h
    Cr = 2000;           % veh/h
    K = 0.04;            % veh/(km*lane)
    v_f = 100;           % km/h
    rho_m = 180;         % veh/(km*lane)
    Dr = 1500;           % veh/h
    a = 1.87;
    lmbda = 4;
    rho_c = 33.5 + E1/3; % veh/(km*lane)

    % Initialize the output matrices
    rho = zeros(4,121);
    v = zeros(4,121);
    wr = zeros(1,121);
    q = zeros(4,1);
    V = zeros(4,1);
    c = zeros(484,1); % Nonlinear inequality constraints
    ceq = []; % Nonlinear equality constraints

    % Initial input
    q0 = 8000 + 100*E1;

    % Initial states
    rho(1:4,1) = 30*ones(4,1);
    v(1:4,1) = 80*ones(4,1);

    for j = 1:1:121
        
        % Change input at k=30
        if (j >= 31)
            q0 = 4000 + 100*E2;
        end
        
        % Calculate traffic flow for each segment
        q(1) = lmbda*(rho(1,j)*v(1,j));
        q(2) = lmbda*(rho(2,j)*v(2,j));
        q(3) = lmbda*(rho(3,j)*v(3,j));
        q(4) = lmbda*(rho(4,j)*v(4,j));
    
        % Calculate the desired speed for each segment
        V(1) = v_f * exp( (-1/a) * ( rho(1,j)/rho_c )^a );
        V(2) = v_f * exp( (-1/a) * ( rho(2,j)/rho_c )^a );
        V(3) = v_f * exp( (-1/a) * ( rho(3,j)/rho_c )^a );
        V(4) = v_f * exp( (-1/a) * ( rho(4,j)/rho_c )^a );

        rj = mapVariable(r(j));
        
        % Calculate the traffic flow that enters segment 4 from the on-ramp
        qr4 = min( [rj*Cr, Dr+(wr(j)/T), Cr*(rho_m-rho(4,j))/(rho_m-rho_c)] );
    
        % Calculate the state update for density
        rho(1,j+1) = rho(1,j) + (T/(lmbda*Li)) * (q0-q(1));
        rho(2,j+1) = rho(2,j) + (T/(lmbda*Li)) * (q(1)-q(2));
        rho(3,j+1) = rho(3,j) + (T/(lmbda*Li)) * (q(2)-q(3));
        rho(4,j+1) = rho(4,j) + (T/(lmbda*Li)) * (q(3)-q(4)+qr4);
        
        % Calculate the state update for velocity
        v(1,j+1) = v(1,j) + (T/tau) * (V(1)-v(1,j)) - (mu*T/(tau*Li)) * (rho(2,j)-rho(1,j)) / (rho(1,j)+K);
        v(2,j+1) = v(2,j) + (T/tau) * (V(2)-v(2,j)) + (T/Li) * v(2,j) * (v(1,j)-v(2,j)) - (mu*T/(tau*Li)) * (rho(3,j)-rho(2,j)) / (rho(2,j)+K);
        v(3,j+1) = v(3,j) + (T/tau) * (V(3)-v(3,j)) + (T/Li) * v(3,j) * (v(2,j)-v(3,j)) - (mu*T/(tau*Li)) * (rho(4,j)-rho(3,j)) / (rho(3,j)+K);
        v(4,j+1) = v(4,j) + (T/tau) * (V(4)-v(4,j)) + (T/Li) * v(4,j) * (v(3,j)-v(4,j));

        % Calculate the state update for queue length
        wr(j+1) = wr(j) + T*(Dr-qr4);

        c((j-1)*4+1,1) = -rho_c + rho(1,j+1);
        c((j-1)*4+2,1) = -rho_c + rho(2,j+1);
        c((j-1)*4+3,1) = -rho_c + rho(3,j+1);
        c((j-1)*4+4,1) = -rho_c + rho(4,j+1);

    end
    
end

% Maps indexes of r to its values
function r = mapVariable(x)
    r_values = .15:.05:.45;
    r = r_values(x);
end

% Cost function
function y = cost_function_ga(r)

    % Group-specific parameters
    E1 = 0 + 1;
    E2 = 7 + 5;

    % Model Parameters
    Li = 1;              % km
    T = 10/3600;         % h
    tau = 18/3600;       % h
    mu = 60;             % km^2/h
    Cr = 2000;           % veh/h
    K = 0.04;            % veh/(km*lane)
    v_f = 100;           % km/h
    rho_m = 180;         % veh/(km*lane)
    Dr = 1500;           % veh/h
    a = 1.87;
    lmbda = 4;
    rho_c = 33.5 + E1/3; % veh/(km*lane)

    % Initialize the output matrices
    rho = zeros(4,122);
    v = zeros(4,122);
    wr = zeros(1,122);
    q = zeros(4,1);
    V = zeros(4,1);

    % Initial input
    q0 = 8000 + 100*E1;

    % Initial states
    rho(1:4,1) = 30*ones(4,1);
    v(1:4,1) = 80*ones(4,1);

    % Initialize output
    y = 0;

    for j = 1:1:121

        % Change input at k=30
        if (j >= 31)
            q0 = 4000 + 100*E2;
        end
        
        % Calculate traffic flow for each segment
        q(1) = lmbda*(rho(1,j)*v(1,j));
        q(2) = lmbda*(rho(2,j)*v(2,j));
        q(3) = lmbda*(rho(3,j)*v(3,j));
        q(4) = lmbda*(rho(4,j)*v(4,j));
    
        % Calculate the desired speed for each segment
        V(1) = v_f * exp( (-1/a) * ( rho(1,j)/rho_c )^a );
        V(2) = v_f * exp( (-1/a) * ( rho(2,j)/rho_c )^a );
        V(3) = v_f * exp( (-1/a) * ( rho(3,j)/rho_c )^a );
        V(4) = v_f * exp( (-1/a) * ( rho(4,j)/rho_c )^a );

        rj = mapVariable(r(j));
        
        % Calculate the traffic flow that enters segment 4 from the on-ramp
        qr4 = min( [rj*Cr, Dr+(wr(j)/T), Cr*(rho_m-rho(4,j))/(rho_m-rho_c)] );
    
        % Calculate the state update for density
        rho(1,j+1) = rho(1,j) + (T/(lmbda*Li)) * (q0-q(1));
        rho(2,j+1) = rho(2,j) + (T/(lmbda*Li)) * (q(1)-q(2));
        rho(3,j+1) = rho(3,j) + (T/(lmbda*Li)) * (q(2)-q(3));
        rho(4,j+1) = rho(4,j) + (T/(lmbda*Li)) * (q(3)-q(4)+qr4);
        
        % Calculate the state update for velocity
        v(1,j+1) = v(1,j) + (T/tau) * (V(1)-v(1,j)) - (mu*T/(tau*Li)) * (rho(2,j)-rho(1,j)) / (rho(1,j)+K);
        v(2,j+1) = v(2,j) + (T/tau) * (V(2)-v(2,j)) + (T/Li) * v(2,j) * (v(1,j)-v(2,j)) - (mu*T/(tau*Li)) * (rho(3,j)-rho(2,j)) / (rho(2,j)+K);
        v(3,j+1) = v(3,j) + (T/tau) * (V(3)-v(3,j)) + (T/Li) * v(3,j) * (v(2,j)-v(3,j)) - (mu*T/(tau*Li)) * (rho(4,j)-rho(3,j)) / (rho(3,j)+K);
        v(4,j+1) = v(4,j) + (T/tau) * (V(4)-v(4,j)) + (T/Li) * v(4,j) * (v(3,j)-v(4,j));
    
        % Calculate the state update for queue length
        wr(j+1) = wr(j) + T * (Dr-qr4);
        
        y = y + T*Li*lmbda*(rho(1,j) + rho(2,j) + rho(3,j) + rho(4,j)) + T*wr(j);

    end

end


% =============================================================
% 7) FUNCTIONS FOR PLOTTING
% =============================================================

% Output optimum states
function [rho, v, wr] = optimum_states(r)
    
    % Group-specific parameters
    E1 = 0 + 1;
    E2 = 7 + 5;
    
    % Model Parameters
    Li = 1;              % km
    T = 10/3600;         % h
    tau = 18/3600;       % h
    mu = 60;             % km^2/h
    Cr = 2000;           % veh/h
    K = 0.04;            % veh/(km*lane)
    v_f = 100;           % km/h
    rho_m = 180;         % veh/(km*lane)
    Dr = 1500;           % veh/h
    a = 1.87;
    lmbda = 4;
    rho_c = 33.5 + E1/3;

    % Initialize the output matrices
    rho = zeros(4,121);
    v = zeros(4,121);
    wr = zeros(1,121);
    q = zeros(4,1);
    V = zeros(4,1);
    
    % Set initial input
    q0 = 8000 + 100*E1;
    
    % Set initial states
    rho(1:4,1) = 30*ones(4,1);
    v(1:4,1) = 80*ones(4,1);
    
    % Calculate optimized values of the variables
    for j = 1:1:120
    
        % Change input at k=30
        if (j >= 31)
            q0 = 4000 + 100*E2;
        end
    
        % Calculate traffic flow for each segment
        q(1) = lmbda*(rho(1,j)*v(1,j));
        q(2) = lmbda*(rho(2,j)*v(2,j));
        q(3) = lmbda*(rho(3,j)*v(3,j));
        q(4) = lmbda*(rho(4,j)*v(4,j));
    
        % Calculate the desired speed for each segment
        V(1) = v_f * exp( (-1/a) * ( rho(1,j)/rho_c )^a );
        V(2) = v_f * exp( (-1/a) * ( rho(2,j)/rho_c )^a );
        V(3) = v_f * exp( (-1/a) * ( rho(3,j)/rho_c )^a );
        V(4) = v_f * exp( (-1/a) * ( rho(4,j)/rho_c )^a );
    
        % Calculate the traffic flow that enters segment 4 from the on-ramp
        qr4 =  r(j)*Cr;
    
        % Calculate the state update for density
        rho(1,j+1) = rho(1,j) + (T/(lmbda*Li)) * (q0-q(1));
        rho(2,j+1) = rho(2,j) + (T/(lmbda*Li)) * (q(1)-q(2));
        rho(3,j+1) = rho(3,j) + (T/(lmbda*Li)) * (q(2)-q(3));
        rho(4,j+1) = rho(4,j) + (T/(lmbda*Li)) * (q(3)-q(4)+qr4);
        
        % Calculate the state update for velocity
        v(1,j+1) = v(1,j) + (T/tau) * (V(1)-v(1,j)) - (mu*T/(tau*Li)) * (rho(2,j)-rho(1,j)) / (rho(1,j)+K);
        v(2,j+1) = v(2,j) + (T/tau) * (V(2)-v(2,j)) + (T/Li) * v(2,j) * (v(1,j)-v(2,j)) - (mu*T/(tau*Li)) * (rho(3,j)-rho(2,j)) / (rho(2,j)+K);
        v(3,j+1) = v(3,j) + (T/tau) * (V(3)-v(3,j)) + (T/Li) * v(3,j) * (v(2,j)-v(3,j)) - (mu*T/(tau*Li)) * (rho(4,j)-rho(3,j)) / (rho(3,j)+K);
        v(4,j+1) = v(4,j) + (T/tau) * (V(4)-v(4,j)) + (T/Li) * v(4,j) * (v(3,j)-v(4,j));
    
        % Calculate the state update for queue length
        wr(j+1) = wr(j) + T*(Dr-qr4);
        
    end

end

% State update function
function [rho_next, v_next, wr_next] = state_update_no_control(x)
    % Model Parameters
    Li = 1;              % km
    T = 10/3600;         % h
    tau = 18/3600;       % h
    mu = 60;             % km^2/h
    Cr = 2000;           % veh/h
    K = 0.04;            % veh/(km*lane)
    v_f = 100;           % km/h
    rho_m = 180;         % veh/(km*lane)
    Dr = 1500;           % veh/h
    a = 1.87;
    lmbda = 4;

    rho = [x(1:4) ; x(4)]; % rho_i for i={1,2,3,4,5} ; rho_4 = rho_5
    v = [x(5); x(5:8)];    % v_i   for i={0,1,2,3,4} ; v_0 = v_1
    wr = x(9);
    rho_c = x(10);
    r = x(11);
    q0 = x(12);
    qr = [0 0 0 min( [r*Cr, Dr+(wr/T), Cr*(rho_m-rho(4))/(rho_m-rho_c)] )];
    q_temp = [q0, 0];
    

    rho_next = zeros(4,1);
    v_next = zeros(4,1);
    
    for i = 1:1:4
        q_temp(2) = lmbda*rho(i)*v(i+1);

        rho_next(i) = rho(i) + ( q_temp(1) - q_temp(2) + qr(i) ) *T/(lmbda*Li);
        
        V = v_f*exp(- ((rho(i)/rho_c)^a) /a);

        v_next(i) = v(i+1) + ( V - v(i+1) )*T/tau + v(i+1)*(v(i) - v(i+1))*T/Li - (mu*T*(rho(i+1) - rho(i))) / (tau*Li* (rho(i)+K) );

        q_temp(1) = q_temp(2);
    end

    wr_next = wr + T*(Dr - qr(4));

end

% =============================================================
% END
% =============================================================