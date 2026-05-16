%% SC42056 Optimization for Systems and Control
% Nonlinear Programming Assignment - Group 22
% Sven Rutgers - 4600150
% Melis Orhun  - 4912071
% =============================================================
% Optimization algorithm with nonlinear equality constraints
% *** BigMac ***
% =============================================================
% *** TABLE OF CONTENTS ***
%  
% 1) PARAMETERS
% 2) OPTIMIZATION
%  * fmincon with SQP and r=0
%  * fmincon with Interior-point and r=0
%  * fmincon with Interior-point and r=0.99
%  * no control case calculation
% 3) PLOTS
%  * TTS, Elapsed Time and Num. of Iterations Table
%  * Plots
% 4) FUNCTIONS
%  * nonlinear constraints function definition
%  * cost function definition
%  * no control case state update function definition
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
E3 = 1 + 0;

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
% 2) OPTIMIZATION
% =============================================================

% Define equality constraints
A_eq = [eye(9) zeros(9,1201)];
b_eq = [30*ones(4,1); 80*ones(4,1); 0];

% Define upper and lower bounds
lb = zeros(1210,1);
ub = [];
ub_i = [rho_c*ones(4,1); inf; inf; inf; inf; inf; 1]; % Set upper bounds of all states for k

% Set initial points
x0 = [];
x0_i = [30*ones(4,1); 80*ones(4,1); 0; 0]; % Set initial points of all states for k

for i = 0:1:120
    ub(i*10+1 : i*10+10,1) = ub_i; % Upper bound matrix for k in [0,120]
    x0(i*10+1 : i*10+10,1) = x0_i; % Set initial points of all states for k in [0,120]
end

% Run different optimization algorithms
% SQP
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); % Use sqp
t_begin_sqp = tic;
[optimum_sqp, out_sqp, exitflag_sqp, data_sqp] = fmincon(@cost_function,x0,[],[],A_eq,b_eq,lb,ub,@nlcon,options);
t_sqp = toc(t_begin_sqp);
% IP
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',inf,'MaxIterations',inf); % Use interior point
t_begin_ip = tic;
[optimum_ip, out_ip, exitflag_ip, data_ip] = fmincon(@cost_function,x0,[],[],A_eq,b_eq,lb,ub,@nlcon,options);
t_ip = toc(t_begin_ip);

% Set initial points
x0 = [];
x0_i = [30*ones(4,1); 80*ones(4,1); 0; 0.99]; % Set initial points of all states for k

for i = 0:1:120
    x0(i*10+1 : i*10+10,1) = x0_i; % Set initial points of all states for k in [0,120]
end

% IP with r=0.99
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',inf,'MaxIterations',inf); % Use interior point
t_begin_ip_099 = tic;
[optimum_ip_099, out_ip_099, exitflag_ip_099, data_ip_099] = fmincon(@cost_function,x0,[],[],A_eq,b_eq,lb,ub,@nlcon,options);
t_ip_099 = toc(t_begin_ip_099);

% Initialize the output matrices
rho_ip = zeros(4,121);
v_ip = zeros(4,121);
wr_ip = zeros(1,121);
r_ip = zeros(1,121);
rho_ip_099 = zeros(4,121);
v_ip_099 = zeros(4,121);
wr_ip_099 = zeros(1,121);
r_ip_099 = zeros(1,121);
rho_sqp = zeros(4,121);
v_sqp = zeros(4,121);
wr_sqp = zeros(1,121);
r_sqp = zeros(1,121);

% Set the optimum points found to the output matrices
for i = 0:1:120
    for j = 1:1:4
        rho_ip(j,i+1) = optimum_ip(i*10 + j);
        v_ip(j,i+1) = optimum_ip(i*10 + j+4);
        rho_ip_099(j,i+1) = optimum_ip_099(i*10 + j);
        v_ip_099(j,i+1) = optimum_ip_099(i*10 + j+4);
        rho_sqp(j,i+1) = optimum_sqp(i*10 + j);
        v_sqp(j,i+1) = optimum_sqp(i*10 + j+4);
    end
    wr_ip(i+1) = optimum_ip(i*10 + 9);
    r_ip(i+1) = optimum_ip(i*10 + 10);
    wr_ip_099(i+1) = optimum_ip_099(i*10 + 9);
    r_ip_099(i+1) = optimum_ip_099(i*10 + 10);
    wr_sqp(i+1) = optimum_sqp(i*10 + 9);
    r_sqp(i+1) = optimum_sqp(i*10 + 10);
end

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
% 3) PLOTS
% =============================================================

% Display data of different algorithms
Tbl = table([out_sqp;out_ip;out_ip_099;y],[t_sqp;t_ip;t_ip_099;0],[data_sqp.iterations;data_ip.iterations;data_ip_099.iterations;0],[exitflag_sqp;exitflag_ip;exitflag_ip_099;0], ...
            'VariableNames',{'Mini TTS (h)','Elapsed Time (s)','Iterations','ExitFlag'},'RowNames',{'BigMac SQP (r=0)','BigMac Interior Point (r=0)','BigMac Interior Point (r=0.99)','No Control Case'});
disp(Tbl)

set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex'); 
set(groot,'defaultLineLineWidth',1)
set(groot,'defaultAxesFontSize',11)
set(groot,'defaultAxesFontWeight',"normal")

% Plot the optimized solutions for BigMac
figure(1) % SQP - density
plot(0:10:1200,rho_sqp(1,1:121),0:10:1200,rho_sqp(2,1:121),0:10:1200,rho_sqp(3,1:121),0:10:1200,rho_sqp(4,1:121))
title(sprintf("Traffic Density vs. Time for Each Lane\n(with SQP and $r_{sp}(k)=0$)"))
xlabel("Time (s)"); ylabel("Traffic Density [veh/(km.lane)]"); grid minor
legend('$\rho_1$(t)','$\rho_2$(t)','$\rho_3$(t)','$\rho_4$(t)')

figure(2) % SQP - mean speed
plot(0:10:1200,v_sqp(1,1:121),0:10:1200,v_sqp(2,1:121),0:10:1200,v_sqp(3,1:121),0:10:1200,v_sqp(4,1:121))
title("Mean Speed vs. Time (with SQP and $r_{sp}(k)=0$)")
xlabel("Time (s)"); ylabel("Mean Speed [km/h]"); grid minor
legend('$v_1$(t)','$v_2$(t)','$v_3$(t)','$v_4$(t)','Location','northwest')

figure(3) % SQP - number of vehicles in queue on the on-ramp
plot(0:10:1200,wr_sqp) 
title(sprintf("Number of Vehicles in Queue on the On-ramp vs. Time\n(with SQP and $r_{sp}(k)=0$)"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor

figure(4) % SQP - ramp metering rate
plot(0:10:1200,r_sqp)
title("Ramp Metering Rate vs. Time (with SQP and $r_{sp}(k)=0$)")
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor

% Plot the optimized solutions for SQP - IP comparison
figure(5) % density comparison
sgtitle("Traffic Density Comparison for SQP and IP Solutions of BigMac","FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,rho_ip(1,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(1,1:121),":");
title("Traffic Density vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{1,ip}$(t)','$\rho_{1,sqp}$(t)')
%
subplot(2,2,2)
plot(0:10:1200,rho_ip(2,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(2,1:121),":");
title("Traffic Density vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{2,ip}$(t)','$\rho_{2,sqp}$(t)')
%
subplot(2,2,3)
plot(0:10:1200,rho_ip(3,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(3,1:121),":");
title("Traffic Density vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{3,ip}$(t)','$\rho_{3,sqp}$(t)')
%
subplot(2,2,4)
plot(0:10:1200,rho_ip(4,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(4,1:121),":");
title("Traffic Density vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$\rho_{4,ip}$(t)','$\rho_{4,sqp}$(t)','Location','southwest')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(6) % speed comparison
sgtitle("Mean Speed Comparison for SQP and IP Solutions of BigMac","FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,v_ip(1,1:121),"--"); hold on; plot(0:10:1200,v_sqp(1,1:121),":");
title("Mean Speed vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{1,ip}$(t)','$v_{1,sqp}$(t)','Location','southeast')
%
subplot(2,2,2)
plot(0:10:1200,v_ip(2,1:121),"--"); hold on; plot(0:10:1200,v_sqp(2,1:121),":");
title("Mean Speed vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{2,ip}$(t)','$v_{2,sqp}$(t)','Location','southeast')
%
subplot(2,2,3)
plot(0:10:1200,v_ip(3,1:121),"--"); hold on; plot(0:10:1200,v_sqp(3,1:121),":");
title("Mean Speed vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{3,ip}$(t)','$v_{3,sqp}$(t)')
%
subplot(2,2,4)
plot(0:10:1200,v_ip(4,1:121),"--"); hold on; plot(0:10:1200,v_sqp(4,1:121),":");
title("Mean Speed vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$v_{4,ip}$(t)','$v_{4,sqp}$(t)')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(7) % queue length comparison
plot(0:10:1200,wr_ip,"--"); hold on; plot(0:10:1200,wr_sqp,":"); 
title("Queue Length vs. Time for SQP and IP Solutions of BigMac")
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor
legend('$w_{r,ip}$(t)','$w_{r,sqp}$(t)')

figure(8) % ramp metering rate comparison
plot(0:10:1200,r_ip,"--"); hold on; plot(0:10:1200,r_sqp,":"); 
title("Ramp Metering Rate vs. Time for SQP and IP Solutions of BigMac")
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor
legend('$r_{ip}$(t)','$r_{sqp}$(t)','Location','southeast')

figure(9) % input q0
subplot(2,1,1)
q_in = [(8000 + 100*E1)*ones(1,30) (4000 + 100*E2)*ones(1,91)];
plot(0:10:1200,q_in)
title("Traffic Flow of Segment 0 ($q_0(t)$) vs. Time")
xlabel("Time (s)"); ylabel("Traffic Flow [veh/h]"); grid minor
%
subplot(2,1,2)
D_r = Dr*ones(1,121);
plot(0:10:1200,D_r)
title("Demand of the On-ramp ($D_r(t)$) vs. Time")
xlabel("Time (s)"); ylabel("Demand [veh/h]"); grid minor

% Plot the optimized and no control solutions
figure(10) % density comparison
sgtitle(sprintf("Traffic Density Comparison\nfor No-Control and SQP-Optimized Solutions of BigMac"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,rho_no_control(1,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(1,1:121),":");
title("Traffic Density vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
%
subplot(2,2,2)
plot(0:10:1200,rho_no_control(2,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(2,1:121),":");
title("Traffic Density vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
%
subplot(2,2,3)
plot(0:10:1200,rho_no_control(3,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(3,1:121),":");
title("Traffic Density vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
%
subplot(2,2,4)
plot(0:10:1200,rho_no_control(4,1:121),"--"); hold on; plot(0:10:1200,rho_sqp(4,1:121),":");
title("Traffic Density vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','northwest')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(11) % speed comparison
sgtitle(sprintf("Mean Speed Comparison\nfor No-Control and SQP-Optimized Solutions of BigMac"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,v_no_control(1,1:121),"--"); hold on; plot(0:10:1200,v_sqp(1,1:121),":");
title("Mean Speed vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','southwest')
%
subplot(2,2,2)
plot(0:10:1200,v_no_control(2,1:121),"--"); hold on; plot(0:10:1200,v_sqp(2,1:121),":");
title("Mean Speed vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','southwest')
%
subplot(2,2,3)
plot(0:10:1200,v_no_control(3,1:121),"--"); hold on; plot(0:10:1200,v_sqp(3,1:121),":");
title("Mean Speed vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','southwest')
%
subplot(2,2,4)
plot(0:10:1200,v_no_control(4,1:121),"--"); hold on; plot(0:10:1200,v_sqp(4,1:121),":");
title("Mean Speed vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$','Location','east')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(12) % queue length comparison
plot(0:10:1200,wr_no_control,"--"); hold on; plot(0:10:1200,wr_sqp,":"); 
title(sprintf("Queue Length vs. Time\nfor for No-Control and SQP-Optimized Solutions of BigMac"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor
legend('$r(k)=1$','$r_{sp}(k)=0$')

% Plot the optimized solutions for IP r={0,0.99}
figure(13) % density comparison
sgtitle(sprintf("Traffic Density Comparison\nfor IP with $r_{sp}(k)=0$ and $r_{sp}(k)=0.99$"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,rho_ip(1,1:121),"--"); hold on; plot(0:10:1200,rho_ip_099(1,1:121),":");
title("Traffic Density vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$')
%
subplot(2,2,2)
plot(0:10:1200,rho_ip(2,1:121),"--"); hold on; plot(0:10:1200,rho_ip_099(2,1:121),":");
title("Traffic Density vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$')
%
subplot(2,2,3)
plot(0:10:1200,rho_ip(3,1:121),"--"); hold on; plot(0:10:1200,rho_ip_099(3,1:121),":");
title("Traffic Density vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$')
%
subplot(2,2,4)
plot(0:10:1200,rho_ip(4,1:121),"--"); hold on; plot(0:10:1200,rho_ip_099(4,1:121),":");
title("Traffic Density vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Traffic Density\n[veh/(km.lane)]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$','Location','southwest')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(15) % speed comparison
sgtitle(sprintf("Mean Speed Comparison\nfor IP with $r_{sp}(k)=0$ and $r_{sp}(k)=0.99$"),"FontSize",12,"FontWeight","bold")
%
subplot(2,2,1)
plot(0:10:1200,v_ip(1,1:121),"--"); hold on; plot(0:10:1200,v_ip_099(1,1:121),":");
title("Mean Speed vs. Time for Lane 1")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$','Location','southeast')
%
subplot(2,2,2)
plot(0:10:1200,v_ip(2,1:121),"--"); hold on; plot(0:10:1200,v_ip_099(2,1:121),":");
title("Mean Speed vs. Time for Lane 2")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$','Location','southeast')
%
subplot(2,2,3)
plot(0:10:1200,v_ip(3,1:121),"--"); hold on; plot(0:10:1200,v_ip_099(3,1:121),":");
title("Mean Speed vs. Time for Lane 3")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$')
%
subplot(2,2,4)
plot(0:10:1200,v_ip(4,1:121),"--"); hold on; plot(0:10:1200,v_ip_099(4,1:121),":");
title("Mean Speed vs. Time for Lane 4")
xlabel("Time (s)"); xlim([0,1200]); ylabel(sprintf("Mean Speed [km/h]")); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$')
set(gcf, 'Position',  [200, 100, 700, 500]);

figure(16) % queue length comparison
plot(0:10:1200,wr_ip,"--"); hold on; plot(0:10:1200,wr_ip_099,":"); 
title(sprintf("Queue Length vs. Time\nfor IP with $r_{sp}(k)=0$ and $r_{sp}(k)=0.99$"))
xlabel("Time (s)"); ylabel("Queue Length [Veh]"); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$')

figure(17) % ramp metering rate comparison
plot(0:10:1200,r_ip,"--"); hold on; plot(0:10:1200,r_ip_099,":"); 
title(sprintf("Ramp Metering Rate vs. Time\nfor IP with $r_{sp}(k)=0$ and $r_{sp}(k)=0.99$"))
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor
legend('$r_{sp}(k)=0$','$r_{sp}(k)=0.99$','Location','southeast')
%
figure(18) % ramp metering rate for no control
plot(0:10:1200,ones(1,121)); 
title(sprintf("Ramp Metering Rate for No Control Case"))
xlabel("Time (s)"); ylabel("Ramp Metering Rate"); ylim("padded"); grid minor

% =============================================================
% 4) FUNCTIONS
% =============================================================

% Nonlinear constraints function
function [c, ceq] = nlcon(x)

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

    % Initialize variables 
    q = zeros(4,1); % Traffic flow
    V = zeros(4,1); % Desired speed for drivers
    c = zeros(120,1); % Nonlinear inequality constraints
    ceq = zeros(1080,1); % Nonlinear equality constraints


    % Initial input
    q0 = 8000 + 100*E1;

    for j = 0:10:1190
        
        % Change input at k=30
        if (j == 300)
            q0 = 4000 + 100*E2;
        end
        
        % Calculate traffic flow for each segment
        q(1) = lmbda*(x(j+1)*x(j+5));
        q(2) = lmbda*(x(j+2)*x(j+6));
        q(3) = lmbda*(x(j+3)*x(j+7));
        q(4) = lmbda*(x(j+4)*x(j+8));
        
        % Calculate the desired speed for each segment
        V(1) = v_f * exp( (-1/a) * ( x(j+1)/rho_c )^a );
        V(2) = v_f * exp( (-1/a) * ( x(j+2)/rho_c )^a );
        V(3) = v_f * exp( (-1/a) * ( x(j+3)/rho_c )^a );
        V(4) = v_f * exp( (-1/a) * ( x(j+4)/rho_c )^a );
    
        % Calculate the traffic flow that enters segment 4 from the on-ramp
        qr4 =  min( [ x(j+10)*Cr, Dr+(x(j+9)/T), Cr*(rho_m-x(j+4))/(rho_m-rho_c) ] );

        % Nonlinear equality constraints on density
        ceq(j+1,1) = -x(j+11) + x(j+1) + (T/(lmbda*Li)) * (q0-q(1));
        ceq(j+2,1) = -x(j+12) + x(j+2) + (T/(lmbda*Li)) * (q(1)-q(2));
        ceq(j+3,1) = -x(j+13) + x(j+3) + (T/(lmbda*Li)) * (q(2)-q(3));
        ceq(j+4,1) = -x(j+14) + x(j+4) + (T/(lmbda*Li)) * (q(3)-q(4)+qr4);

        % Nonlinear equality constraints on velocity
        ceq(j+5,1) = -x(j+15) + x(j+5) + (T/tau) * (V(1)-x(j+5)) - (mu*T/(tau*Li)) * (x(j+2)-x(j+1)) / (x(j+1)+K);
        ceq(j+6,1) = -x(j+16) + x(j+6) + (T/tau) * (V(2)-x(j+6)) + (T/Li) * x(j+6) * (x(j+5)-x(j+6)) - (mu*T/(tau*Li)) * (x(j+3)-x(j+2)) / (x(j+2)+K);
        ceq(j+7,1) = -x(j+17) + x(j+7) + (T/tau) * (V(3)-x(j+7)) + (T/Li) * x(j+7) * (x(j+6)-x(j+7)) - (mu*T/(tau*Li)) * (x(j+4)-x(j+3)) / (x(j+3)+K);
        ceq(j+8,1) = -x(j+18) + x(j+8) + (T/tau) * (V(4)-x(j+8)) + (T/Li) * x(j+8) * (x(j+7)-x(j+8));
    
        % Nonlinear equality constraint on queue length
        ceq(j+9,1) = -x(j+19) + x(j+9) + T*(Dr-qr4);
        
        %c(j+1,1) = -0.005 + ( x(j+20) - x(j+10) )^2; % Extra term for smooth graphs
    
    end
    
end

% Cost function
function y = cost_function(x)
    Li = 1;              % km
    T = 10/3600;         % h
    lmbda = 4;
    y = 0;

    for i = 0:10:1200
        y = y + T*Li*lmbda*(x(i+1) + x(i+2) + x(i+3) + x(i+4)) + T*x(i+9);% + 0*( x(i+10)-x(i-1) )^2; % Extra term for smooth graphs
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