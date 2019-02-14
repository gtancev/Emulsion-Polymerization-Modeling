function [ dm ] = emuslion1( t,n,n_s,n_M_0,V_p_0 )

k_p = 715; % L/(mol*s)
n_bar = 0.5;
N_A = 6.022*10^23; % 1/mol
M_M = 0.1; % kg/mol
p_M = 0.94; % kg/L
p_P = 1.1; % kg/L

%%

n_M = n(1); % amount of total monomer available
n_M_p = n(2); % amount of monomer in particles
n_M_d = n_M-n_M_p; % amount of monomer in dispersed phase

%%

v_pol = V_p_0+(n_M_0-n_M)*M_M/(n_s*p_P);
v_p = v_pol+(n_M_p/n_s)*(M_M/p_M);
phi = (v_p-v_pol)/v_p;
c_M = n_M_p/(n_s*v_p);

%%

dn_M = - k_p*c_M*n_bar/N_A*n_s;
dv_pol = - dn_M*M_M/(p_P*n_s);

%%
    
    if n_M_d > 0
        
        dn_M_P = n_s.*c_M.*1./(1-phi).*dv_pol;
    
    else

        dn_M_P = dn_M;
        
    end
    
%%

dm = [dn_M; dn_M_P];

end

