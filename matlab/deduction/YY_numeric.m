clear; close all; clc;

load('f_YY_menor60.mat')
load('f_YY_harmonic_menor60.mat')
load('f_YY_maior60.mat')
load('f_YY_harmonic_maior60.mat')



addpath('utils');
addpath('data');


% wire data
awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32
awg_chosen = 38;
strand = 180;
d_l = awg_wires(awg_chosen)*1e-3; % [m]



Vp = 400;
L = 61e-6;
n = 5/9;
d = 1;
fs = 100e3;
phi = 67*pi/180;
k = 1;


% parametros trafo e indutor (numero de voltas e espiras em paralelo)
NL = 14;
Np = 9;
Ns = 5;
Sp = 2;
Ss = 3;

% Vp = double(Vp_num);
% L = double(L_num);
% n = double(n_num);
% d = double(d_num);
% fs = double(fs_num);
% phi = double(phi_num);


if (phi<pi/3 && phi>0)
    M = f_YY_menor60(Vp,L,n,d,fs,phi,1);
else
    M = f_YY_maior60(Vp,L,n,d,fs,phi,1);
end

I_in_med = M(1);
I_in_rms = M(2);
I_out_med = M(3);
I_out_rms = M(4);
I_L_rms = M(5);
I_trf_rms = M(6);
I_sw_p_rms = M(7);
I_sw_s_rms = M(8);
I_sw_p_on = M(9);
I_sw_s_on = M(10);
Wb_trf = M(11);
Wb_indutor = M(12);
P_med = I_in_med*Vp;

%% Conduction loss

Pcond = f_cond_loss(I_sw_p_rms,I_sw_s_rms);


[ZVS_p, ZVS_s, p_switch, s_switch] = f_switch_loss(I_sw_p_on,I_sw_s_on,Vp,d,fs);

%% Copper losses
Pind_cu = 0;
Ptrf_cu = 0;
k_max = 5;
for i=1:20
%     k = i
    if (phi<pi/3)
        M_harmonic = f_YY_harmonic_menor60(Vp,L,n,d,fs,phi,i)/sqrt(2);
    else
        M_harmonic = f_YY_harmonic_maior60(Vp,L,n,d,fs,phi,i)/sqrt(2);
    end

    ik_L = M_harmonic(1);
    ik_trf = M_harmonic(2);

    Rk_ind = f_get_resistance_ind(fs*i,NL,d_l);
    Rk_trf = f_get_resistance_trf(fs*i,Np,Ns,Sp,Ss,d_l);

    Pind_cu = Pind_cu + Rk_ind*ik_L^2;
    Ptrf_cu = Ptrf_cu + Rk_trf*ik_trf^2;
end

Pcu = 3*Pind_cu + Ptrf_cu; % Copper losses


%% perdas no nucleo do indutor


omega = 2*pi*fs;
a_L = 1.394;
b_L = 2.248;
kc_L = 2.448;
int_ki = integral(@(theta) abs(cos(theta)).^a_L,0,2*pi);
ki_L = kc_L/(2^(b_L-a_L)*(2*pi)^(a_L-1)*int_ki);

AeL = 421.3e-6;
Ve_L = 52.1e-6;


if (phi<pi/3 && phi > 0)
    integrais_L = (3 ^ (-a_L)) * NL ^ (-a_L) * Vp ^ a_L * AeL ^ (-a_L) * ((2 ^ a_L + 2) * (pi - (3 * phi)) * (abs(-n + d) ^ a_L) / 0.3e1 + (phi * (abs(-n + 2 * d) ^ a_L + (n + d) ^ a_L + abs(-2 * n + d) ^ a_L))) * (n ^ (-a_L)) / omega;
else
    integrais_L = 0.2e1 / 0.3e1 * (3 ^ (-a_L)) * AeL ^ (-a_L) * NL ^ (-a_L) * ((pi - 0.3e1 / 0.2e1 * phi) * (abs(-2 * n + d) ^ a_L) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * (abs(-n + d) ^ a_L) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * ((2 * n + d) ^ a_L) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * ((n + 2 * d) ^ a_L) + (pi - 0.3e1 / 0.2e1 * phi) * (abs(-n + 2 * d) ^ a_L + (n + d) ^ a_L)) * (n ^ (-a_L)) * Vp ^ a_L / omega;
end


deltaB_L = Wb_indutor/(NL*AeL);
T = 1/fs;
PvL = ki_L*2/T*abs(deltaB_L)^(b_L-a_L)*integrais_L;

Pc_L = 3*Ve_L*PvL;

%% perdas no nucleo do transformador

a_trf = 1.585;
b_trf = 2.756;
kc_trf = 0.402;
int_ki = integral(@(theta) abs(cos(theta)).^a_trf,0,2*pi);
ki_trf = kc_trf/(2^(b_trf-a_trf)*(2*pi)^(a_trf-1)*int_ki);

Ae_trf = 683e-6;
Ve_trf = 146.79e-6;


integrais_trf = (3 ^ (-1 - a_trf)) / omega * Vp ^ a_trf * d ^ a_trf * n ^ (-a_trf) * Ns ^ (-a_trf) * Ae_trf ^ (-a_trf) * pi * (2 ^ a_trf + 2);

deltaB_trf = Wb_trf/(Ns*Ae_trf);
Pv_trf = ki_trf*2/T*abs(deltaB_trf)^(b_trf-a_trf)*integrais_trf;

Pc_trf = Ve_trf*Pv_trf


%% calculo eficiencia

Ptotal = Pcond + Pcu + p_switch + s_switch + Pc_L + Pc_trf;

eficiency = (P_med - Ptotal)/P_med

