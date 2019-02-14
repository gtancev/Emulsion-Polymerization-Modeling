function [ eqs ] = emulsion3(t, x, N_M_zero, nr_part, k_d, k_p, k_e, k_t, k_des, phi_M_star, n_mean, Avogadro, MW_M, rho_M, rho_P, rho_S, v_pol_zero)

N_M = x(1);
N_MP = x(2);

N_MD = N_M  - N_MP;

% polymer and particle volume 
v_pol = v_pol_zero + (N_M_zero - N_M).*MW_M./(nr_part.*rho_P);
v_p = v_pol + N_MP./nr_part.*MW_M./rho_M;
v_m = v_p - v_pol;
phi = (v_p-v_pol)./v_p;

r_part = (6.*v_p/pi).^(1/3)/2/10; % meter

% conc. of M in the particle
M_P = N_MP./(nr_part.*v_p);

d_N_M = -  k_p.*M_P.*n_mean(r_part)./Avogadro.*nr_part;

d_v_pol = - d_N_M.*MW_M/(rho_P.*nr_part);

% distinguish phases
    % growth
    if N_MD > 0
        d_N_MP = nr_part.*M_P.*1./(1-phi).*d_v_pol;
%          = N_p * M_p * (-1/N_p * MW_M / rho_p) * dy(1) * 1/(1-phi); % dN_p

    % monomer depletion    
    else

        d_N_MP = d_N_M;
        
    end

eqs = [d_N_M;d_N_MP;];
% keyboard
end

