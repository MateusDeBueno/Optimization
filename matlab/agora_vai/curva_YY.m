clear; close all; clc;

addpath('utils')
addpath('utils_transf')
addpath('utils_loss')
addpath('data');

awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32

l.eq = load('YY.mat');
l.sw = load('f_fitted_off.mat');
l.sw = load('f_fitted_on.mat');
l.sw.Ronp = 90e-3;
l.sw.Rons = 90e-3;


l.L.a = 1.394;
l.L.b = 2.248;
l.L.kc = 2.448;
l.L.int_ki = integral(@(theta) abs(cos(theta)).^l.L.a,0,2*pi);
l.L.ki = l.L.kc/(2^(l.L.b-l.L.a)*(2*pi)^(l.L.a-1)*l.L.int_ki);
l.L.Ac = 421.3e-6;
l.L.Ve = 52.1e-6;
l.L.N = 14;
l.L.strand = 180;
l.L.awg = 38;
l.L.dl = awg_wires(l.L.awg)*1e-3; % [m]
l.L.MLT = 140*1e-3; % [m]
l.L.Nl = 1;

l.tr.a = 1.585;
l.tr.b = 2.756;
l.tr.kc = 0.402;
l.tr.int_ki = integral(@(theta) abs(cos(theta)).^l.tr.a,0,2*pi);
l.tr.ki = l.tr.kc/(2^(l.tr.b-l.tr.a)*(2*pi)^(l.tr.a-1)*l.tr.int_ki);
l.tr.Ac = 683e-6;
l.tr.Ve = 146.79e-6;
l.tr.Np = 9;
l.tr.Ns = 5;
l.tr.Sp = 2;
l.tr.Ss = 3;
l.tr.N = l.tr.Ns;
l.tr.strand = 180;
l.tr.awg = 38;
l.tr.dl = awg_wires(l.tr.awg)*1e-3; % [m]
l.tr.MLT = 230*1e-3; % [m]

l.sC.Rac100 = 4.3e-3;


l.pr.dt = 440e-9;
% l.pr.phi = deg2rad(67);
l.pr.Vi = 400;
l.pr.d = 1;
l.pr.fs = 100e3;
l.pr.Ldab = 61e-6;
l.pr.Ld1 = 1e-12;
l.pr.n = l.tr.Ns/l.tr.Np;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 700000000e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);





ang_min = 40*pi/180;
ang_max = 70*pi/180;
ii = 0;
for phi=ang_min:(ang_max-ang_min)/30:ang_max
    ii=ii+1;
    l.pr.phi = phi;
    l.pr.phi
    vector_phi(ii) = phi*180/pi;
    [efc(ii),ptot(ii),Pm(ii),cSw_p(ii),cSw_s(ii),sSw_p(ii),sSw_s(ii),Pil_cu(ii),PiL_cu(ii),PiLd_cu(ii),Psc(ii),P_core_L(ii),P_core_tr(ii)] = YY_loss(l);
    l.pr.phi
end

figure
scatter(Pm,efc,'filled')
grid on

phi_exp = [67 64 60 50 40];
n_exp = [94.56 94.63 94.23 93.78 92.62];

figure
hold on
scatter(vector_phi,efc,'filled')
scatter(phi_exp,n_exp/100,'filled')
hold off
grid on


figure
scatter(vector_phi,Pm,'filled')
grid on

sw_p = cSw_p+sSw_p;
sw_s = cSw_s+sSw_s;
pcu_tr = Pil_cu+PiL_cu;
pcu_L = PiLd_cu;
P_core_L;
P_core_tr;
figure
bar(vector_phi,[sw_p',sw_s',pcu_tr',pcu_L',P_core_L',P_core_tr'])
grid on
legend({'sw_p','sw_s','pcu_tr','pcu_L','P_core_L','P_core_tr'})

