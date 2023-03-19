clear; close all; clc;

addpath('utils')
addpath('utils_transf')
addpath('utils_loss')
addpath('data');

awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32

l.eq = load('YY.mat');
l.eq = load('DiY.mat');
l.sw = load('f_fitted_off.mat');
l.sw = load('f_fitted_on.mat');
l.sw.Ronp = 102e-3;
l.sw.Rons = 102e-3;


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
l.tr.kc = 0.522;
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
l.C.Rb_ac10k = 4.2e-3;


l.pr.dt = 0e-9;
l.pr.phi = deg2rad(-27);
l.pr.Vi = 400;
l.pr.d = 1;
l.pr.fs = 100e3;
l.pr.Ldab = 61e-6;
l.pr.Ld1 = 1.4e-6;
l.pr.n = l.tr.Ns/l.tr.Np;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 700e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);

%% equations
% out = f_equationsYY(l);
out = f_equationsDiY(l);
C = num2cell(out);
[hbrm,HBrm,Ip,Is,iiRMS,iiME,ioRMS,ioME,Pm,idrm,ilrm,iLrm,iSwPrm,iSwSrm,P_core_tr,Bpk_tr,P_core_L,Bpk_L] = C{:};

%% semi conductor loss

cSw_p = 6*l.sw.Ronp*iSwPrm^2;
cSw_s = 6*l.sw.Rons*iSwSrm^2;

sSw_p = 6*f_sw_loss(l.pr.fs, l.pr.Vi, Ip, l);
sSw_s = 6*f_sw_loss(l.pr.fs, l.pr.Vi*l.pr.d, Is, l);

%% wire loss
Pil_cu = 0;
PiL_cu = 0;
PiLd_cu = 0;
for nn=1:100
    [R_ac_p, R_ac_s] = f_get_resistance_trf(nn,0.8,l);
    RacL = f_get_resistance_ind(nn,0.8,l);
    
    out = fh_equationsDiY(l,nn);
    C = num2cell(out);
    [ilrm_cn,iLrm_cn,idrm_cn] = C{:};
    
    Pil_cu = Pil_cu + R_ac_p*ilrm_cn^2;
    PiL_cu = PiL_cu + R_ac_s*iLrm_cn^2;
    PiLd_cu = PiLd_cu + RacL*idrm_cn^2;
end

Pil_cu = 3*Pil_cu;
PiL_cu = 3*PiL_cu;
PiLd_cu = 3*PiLd_cu;

%% serie capacitor
Psc = 3*l.sC.Rac100*hbrm^2 + 3*l.sC.Rac100*HBrm^2;

%% bus capacitor
ic_i = sqrt(iiRMS^2 - iiME^2);
Pbusi = ic_i*ic_i*l.C.Rb_ac10k;

ic_o = sqrt(ioRMS^2 - ioME^2);
Pbuso = ic_o*ic_o*l.C.Rb_ac10k;
%%

ptot = cSw_p+cSw_s+sSw_p+sSw_s+Pil_cu+PiL_cu+PiLd_cu+3*P_core_L+P_core_tr+Psc + Pbusi + Pbuso;

Pm
efc = (Pm-ptot)/Pm


