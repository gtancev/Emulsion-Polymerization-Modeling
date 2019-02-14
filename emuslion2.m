function [ dm ] = emuslion2( t,n,n_s,n_M_0,V_p_0,V )

k_p = 715; % L/(mol*s)
k_des = 1e-1; % 1/s
k_e = 1e-6; % dm/s
k_t = 9.8e6; % L/(mol*s)
k_d = 5.55e-6; % 1/s

N_A = 6.022*10^23; % 1/mol
M_M = 0.1; % kg/mol
p_M = 0.94; % kg/L
p_P = 1.10; % kg/L
w_I = 0.01;
m_M_0 = n_M_0*M_M; % kg
m_I = w_I*m_M_0; % kg
M_I = 0.164; % kg/mol
f = 0.5;
c_I = m_I/(M_I*V); % mol/L
c_R = sqrt(2*f*k_d*c_I/k_t); % mol/L

%%

n_M = n(1); % amount of total monomer available
n_M_p = n(2); % amount of monomer in particles
n_M_d = n_M-n_M_p; % amount of monomer in dispersed phase

%%

v_pol = V_p_0+(n_M_0-n_M)*M_M/(n_s*p_P); 
v_p = v_pol+(n_M_p/n_s)*(M_M/p_M);
phi = (v_p-v_pol)/v_p;
c_M = n_M_p/(n_s*v_p);

diam = (6*v_p/pi)^(1/3);
A_p = 4*pi*(diam/2)^2;
rho = A_p*k_e*c_R*N_A;
n_bar = (rho/(2*rho+k_des));

%%

dn_M = - k_p*c_M*n_bar/N_A*n_s;
dv_pol = - dn_M*M_M/(p_P*n_s);

%%
    
    if n_M_d > 0
        
        dn_M_P = n_s*c_M*1/(1-phi)*dv_pol;
    
    else

        dn_M_P = dn_M;
        
    end
    
%%

dm = [dn_M; dn_M_P];

end