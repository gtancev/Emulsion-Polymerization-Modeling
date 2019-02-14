%% Polymerization Reaction & Colloid Engineering
% George Tancev

clear all; close all; clc;
options = optimset('Display','none','MaxFunEvals',1000);

%% data

V = 10; % L
d_p_0 = 50e-8; % dm
w_s = 0.05;
c_M_0 = 2; % mol/L
n_M_0 = c_M_0*V;
w_I_0 = 0.01;
a_E = 50*(10^(-9))^2; % dm^2
N_A = 6.022*10^23; % 1/mol

M_M = 0.1; % kg/mol
M_I = 0.164; % kg/mol
p_M = 0.94; % kg/L
p_P = 1.1; % kg/L
p_W = 1.0; % kg/L
c_M = p_M/M_M;
phi_M = 0.5;
w_I = 0.01;
f = 0.5;

k_p = 715; % L/(mol*s)
k_des = 1e-1; % 1/s
k_e = 1e-6; % dm/s
k_t = 9.8*10^6; % L/(mol*s)
k_d = 5.55e-6; % 1/s

cmc = 40e-3; % mol/L
K_E = 100; % L/mol

m_M_0 = n_M_0*M_M; % initial mass of monomer (kg)
V_M_0 = m_M_0/p_M; % initial volume of monomer (L)

%% a)

mass_water = @(m_W_0)(V-p_M*m_M_0-p_W*m_W_0-p_P*(w_s*(m_M_0+m_W_0)/(1-w_s)));

m_W_0 = fsolve(mass_water,0.5,options); % mass of water (kg)
V_W_0 = m_W_0/p_W; % initial volume of water (L)
m_tot = m_M_0+m_W_0; % total mass of solution (kg) without seeds
m_s = w_s*(m_tot)/(1-w_s); % mass of seeds (kg)

V_p_0 = (4/3)*pi*(d_p_0/2)^3; % volume of a seed
m_1_s = V_p_0*p_P; % mass of one seed without monomer

n_s = m_s/m_1_s; % number of seeds
V_s_tot = n_s*2*V_p_0; % total volume of swollen seed
d_p_m = (12*V_p_0/pi)^(1/3); % diameter of swollen seed
A_s_tot = n_s*4*pi*(d_p_m/2)^2; % total surface of all particles at t = 0
A_1_s = 4*pi*(d_p_m/2)^2; % surface of one particle at t = 0
V_1_s = 2*V_p_0; % volume of one particle at t = 0

F_E_inf = (A_1_s/V_1_s)/(a_E*6.02*10^23); % maximum amount of emulsifier adsorbed in mol/vol

F_E_min = 0.2*F_E_inf; % min. adsorbed amount per particle
E_w_min = (F_E_min/F_E_inf)/(K_E*(1-F_E_min/F_E_inf));

E_w_max = cmc;
F_E_max = F_E_inf*(K_E*E_w_max)/(1+K_E*E_w_max); % max. adsorbed amount per particle

n_0_min = E_w_min*V_W_0+F_E_min*n_s*V_1_s; % minimum amount in mol
n_0_max = E_w_max*V_W_0+F_E_max*n_s*V_1_s; % maximum amount in mol

m_I = w_I*m_M_0;
c_I = m_I/(M_I*V);
c_R = sqrt(2*f*k_d*c_I/k_t);

%% b)

m0 = [n_M_0 n_s*V_p_0*c_M];
tspan = 0:1:3600;
options2 = odeset('RelTol',1e-6,'AbsTol',1e-10);
[t,m] = ode15s(@(t,n)emulsion1( t,n,n_s,n_M_0,V_p_0 ),tspan,m0,options2);

v_pol = V_p_0+(n_M_0-m(:,1))./n_s*(M_M/p_P);
v_p = v_pol+(m(:,2)./n_s)*(M_M/p_M);

phi = (v_p-v_pol)./v_p; % phi
diam = (6.*v_p./pi).^(1/3); % diameter over time (nm)
V_pol = n_s*v_pol; % total volume of polymer
V_p = n_s*v_p;
conv1 = (n_M_0-m(:,1))/n_M_0;
A_p = 4.*pi.*(diam./2).^2;
rho = A_p.*k_e.*c_R.*N_A; % 1/s
n_bar = (rho./(2.*rho));

c_M_1 = m(:,2)./(n_s.*v_p);
c_R_1 = n_bar./(N_A.*v_p);
tau_p1 = 1./(k_p.*c_M_1);
tau_in1 = 1./rho;
n_N = (tau_in1./tau_p1);
n_W = 2.*n_N;

figure(1);
subplot(6,1,1);
plot(t./3600,conv1);
title('conversion vs. time');
xlabel('time / [h]');
ylabel('conversion');

subplot(6,1,2);
plot(conv1,V_pol,conv1,V_p);
axis([0 1 0 2.5]);
title('total volume of polymer and particles vs. conversion');
xlabel('x_M');
ylabel('V / [L]');
legend('V_{pol}','V_p','Location','best');

subplot(6,1,3);
plot(conv1,phi);
axis([0 1 0 0.55]);
title('\phi vs. conversion');
xlabel('x_M');
ylabel('\phi');

subplot(6,1,4);
plot(conv1,diam.*10^8);
axis([0 1 50 100]);
title('diameter of a particle vs. conversion');
xlabel('x_M');
ylabel('d / [nm]');

subplot(6,1,5);
plot(t/3600,n_bar);
axis([0 1 0 1]);
title('n vs. conversion');
xlabel('x_M');
ylabel('n');

subplot(6,1,6);
plot(conv1,n_N,conv1,n_W);
title('number and weight average vs. conversion');
xlabel('x_M');
ylabel('n_N, n_W');
legend('n_N','n_W','Location','best');

%% c)

m0 = [n_M_0 n_s*V_p_0*c_M];
tspan2 = 0:1:3600;
[t2,m2] = ode15s(@(t,n)emulsion2( t,n,n_s,n_M_0,V_p_0,V ),tspan2,m0,options2);

v_pol2 = V_p_0+(n_M_0-m2(:,1))./n_s*(M_M/p_P);
v_p2 = v_pol2+(m2(:,2)./n_s)*(M_M/p_M);

phi = (v_p2-v_pol2)./v_p2; % phi
diam2 = (6*v_p2./pi).^(1/3); % diameter over time (dm)
V_pol = n_s*v_pol2; % total volume of polymer
V_p = n_s*v_p2;
conv2 = (n_M_0-m2(:,1))/n_M_0;
A_p2 = 4.*pi.*(diam2./2).^2;
rho = A_p2.*k_e.*c_R.*N_A; % 1/s
n_bar = (rho./(2.*rho+k_des));

c_M_2 = m2(:,2)./(n_s.*v_p2);
c_R_2 = n_bar./(N_A.*v_p2);
tau_p = 1./(k_p.*c_M_2);
tau_in = 1./rho;
tau_out = 1./(k_des.*n_bar);

n_N = 1./(tau_p./tau_in+tau_p./tau_out);
n_W = 2.*n_N;

figure(2);
subplot(6,1,1);
plot(t2./3600,conv2);
title('conversion vs. time');
xlabel('time / [h]');
ylabel('conversion');

subplot(6,1,2);
plot(conv2,V_pol,conv2,V_p);
axis([0 1 0 2.5]);
title('total volume of polymer and particles vs. conversion');
xlabel('x_M');
ylabel('V / [L]');
legend('V_{pol}','V_p','Location','best');

subplot(6,1,3);
plot(conv2,phi);
axis([0 1 0 0.55]);
title('\phi vs. conversion');
xlabel('x_M');
ylabel('\phi');

subplot(6,1,4);
plot(conv2,diam2.*10^8);
axis([0 1 50 100]);
title('diameter of a particle vs. conversion');
xlabel('x_M');
ylabel('d / [nm]');

subplot(6,1,5);
plot(t2/3600,n_bar);
axis([0 1 0 0.5]);
title('n vs. conversion');
xlabel('time / [h]');
ylabel('n');

subplot(6,1,6);
plot(conv2,n_N,conv2,n_W);
title('number and weight average vs. conversion');
xlabel('x_M');
ylabel('n_N, n_W');
legend('n_N','n_W','Location','best');

%%
