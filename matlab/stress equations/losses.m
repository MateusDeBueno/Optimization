clear; close all; clc;

addpath('utils')
addpath('utils_transf')

load('trafo.mat')

syms phi real
syms fs Vi d Ldab Ld1 Ld2 Lm n Po dt t real positive


Vi_num = 400;
d_num = 1;
Vo_num = Vi_num*d_num;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 10e-6;
Ld2_num = 6e-6;
Lm_num = 700e-6;
phi_num = deg2rad(50);  %[MUDAR]
n_num = 5/9;
Po_num = 2000;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
% k_num = M_num/sqrt(L1_num*L2_num);


%% gerar graficos

out = trafo.DinD.str.f3;

out(Ldab_num,n_num,Ld1_num,Ld2_num,Lm_num,phi_num,fs_num,Vi_num,Vo_num)

% [IL_rms,Itrfsec_rms,Iin_med,Iin_rms,Iout_med,Iout_rms,Isw_p_rms,Isw_s_rms,Ip,Is]
% 
% out()
%%

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
Ld2_num = 2e-6;
Lm_num = 700e-6;
% phi_num = deg2rad(67);  %[MUDAR]
n_num = 5/9;
Vo_num = d_num*Vi_num;
k_num = 1;
P_num = 3000;

vec_phi = -180:360/1000000:180;
Iout_med = zeros(1,1001);
count = 0;
for phi_num=vec_phi
    count = count + 1;
    [out] = fYY(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, deg2rad(phi_num), fs_num, Vi_num, Vo_num);
    Pout_med(count) = out(5)*Vo_num;
end

[k] = get_k(Pout_med,P_num);
vec_phi(k)
phi_num = deg2rad(vec_phi(k));
%%



%%

% % % % 
% % % % %%
% % % % [out] = fYY(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num);
% % % % 
% % % % IL_rms = out(1)
% % % % Itrfsec_rms = out(2)
% % % % Iin_med = out(3)
% % % % Iin_rms = out(4)
% % % % Iout_med = out(5)
% % % % Iout_rms = out(6)
% % % % Isw_p_rms = out(7)
% % % % Isw_s_rms = out(8)
% % % % Ip = out(9)
% % % % Is  = out(10)
% % % % 
% % % % 
% % % % 
% % % % %%
% % % % 
% % % % 
% % % % 
% % % % % wire data
% % % % awg_data = readtable('awg_table.txt');
% % % % awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32
% % % % awg_chosen = 38;
% % % % strand = 180;
% % % % d_l = awg_wires(awg_chosen)*1e-3; % [m]
% % % % 
% % % % NL = 14;
% % % % Np = 9;
% % % % Ns = 5;
% % % % Sp = 2;
% % % % Ss = 3;
% % % % 
% % % % Pcond = f_cond_loss(Isw_p_rms,Isw_s_rms)
% % % % 
% % % % 
% % % % [ZVS_p, ZVS_s, p_switch, s_switch] = f_switch_loss(I_sw_p_on,I_sw_s_on,Vp,d,fs);
% % % %     
% % % % %% Copper losses
% % % % Pind_cu = 0;
% % % % Ptrf_cu = 0;
% % % % k_max = 5;
% % % % for i=1:50
% % % % %     k = i
% % % %     if (phi<pi/3)
% % % %         M_harmonic = f_YY_harmonic_menor60(Vp,L,n,d,fs,phi,i)/sqrt(2);
% % % %     else
% % % %         M_harmonic = f_YY_harmonic_maior60(Vp,L,n,d,fs,phi,i)/sqrt(2);
% % % %     end
% % % %     
% % % %     ik_L = M_harmonic(1);
% % % %     ik_trf = M_harmonic(2);
% % % %     
% % % %     Rk_ind = f_get_resistance_ind(fs*i,NL,d_l);
% % % %     Rk_trf = f_get_resistance_trf(fs*i,Np,Ns,Sp,Ss,d_l);
% % % %     
% % % %     Pind_cu = Pind_cu + Rk_ind*ik_L^2;
% % % %     Ptrf_cu = Ptrf_cu + Rk_trf*ik_trf^2;
% % % % end
% % % %     
% % % % Pcu = 3*Pind_cu + Ptrf_cu; % Copper losses
% % % %     
% % % %     
% % % % %% perdas no nucleo do indutor
% % % % 
% % % % 
% % % % % omega = 2*pi*fs;
% % % % % a_L = 1.394;
% % % % % b_L = 2.248;
% % % % % kc_L = 2.448;
% % % % % int_ki = integral(@(theta) abs(cos(theta)).^a_L,0,2*pi);
% % % % % ki_L = kc_L/(2^(b_L-a_L)*(2*pi)^(a_L-1)*int_ki);
% % % % % 
% % % % % AeL = 421.3e-6;
% % % % % Ve_L = 52.1e-6;
% % % % % 
% % % % % 
% % % % % if (phi<pi/3 && phi > 0)
% % % % %     integrais_L = (3 ^ (-a_L)) * NL ^ (-a_L) * Vp ^ a_L * AeL ^ (-a_L) * ((2 ^ a_L + 2) * (pi - (3 * phi)) * (abs(-n + d) ^ a_L) / 0.3e1 + (phi * (abs(-n + 2 * d) ^ a_L + (n + d) ^ a_L + abs(-2 * n + d) ^ a_L))) * (n ^ (-a_L)) / omega;
% % % % % else
% % % % %     integrais_L = 0.2e1 / 0.3e1 * (3 ^ (-a_L)) * AeL ^ (-a_L) * NL ^ (-a_L) * ((pi - 0.3e1 / 0.2e1 * phi) * (abs(-2 * n + d) ^ a_L) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * (abs(-n + d) ^ a_L) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * ((2 * n + d) ^ a_L) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * ((n + 2 * d) ^ a_L) + (pi - 0.3e1 / 0.2e1 * phi) * (abs(-n + 2 * d) ^ a_L + (n + d) ^ a_L)) * (n ^ (-a_L)) * Vp ^ a_L / omega;
% % % % % end
% % % % %    
% % % % %     
% % % % % deltaB_L = Wb_indutor/(NL*AeL);
% % % % % T = 1/fs;
% % % % % PvL = ki_L*2/T*abs(deltaB_L)^(b_L-a_L)*integrais_L;
% % % % %     
% % % % % Pc_L = 3*Ve_L*PvL;
% % % %     
% % % % %% perdas no nucleo do transformador
% % % % 
% % % % % a_trf = 1.585;
% % % % % b_trf = 2.756;
% % % % % kc_trf = 0.402;
% % % % % int_ki = integral(@(theta) abs(cos(theta)).^a_trf,0,2*pi);
% % % % % ki_trf = kc_trf/(2^(b_trf-a_trf)*(2*pi)^(a_trf-1)*int_ki);
% % % % % 
% % % % % Ae_trf = 683e-6;
% % % % % Ve_trf = 146.79e-6;
% % % % % 
% % % % % 
% % % % % integrais_trf = (3 ^ (-1 - a_trf)) / omega * Vp ^ a_trf * d ^ a_trf * n ^ (-a_trf) * Ns ^ (-a_trf) * Ae_trf ^ (-a_trf) * pi * (2 ^ a_trf + 2);
% % % % % 
% % % % % deltaB_trf = Wb_trf/(Ns*Ae_trf);
% % % % % Pv_trf = ki_trf*2/T*abs(deltaB_trf)^(b_trf-a_trf)*integrais_trf;
% % % % % 
% % % % % Pc_trf = Ve_trf*Pv_trf;
% % % % 
% % % % 
% % % % %% calculo eficiencia
% % % %     
% % % % Ptotal = Pcond + Pcu + p_switch + s_switch + Pc_L + Pc_trf;
% % % %    
% % % % eficiency = (P_med - Ptotal)/P_med;