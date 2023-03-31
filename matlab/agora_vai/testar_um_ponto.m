clear; close all; clc;

% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


addpath('utils')
addpath('utils_transf')
addpath('utils_loss')
addpath('data');
addpath('histcompress')
addpath('consolidator')

awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32

load('YY.mat');
load('YD.mat');
load('DfD.mat');
load('DiD.mat');
load('DiY.mat');
load('DfY.mat');

% Combine the structs into a single struct
l.eq = struct('YY', YY, 'YD', YD, 'DfD', DfD, ...
              'DiD', DiD, 'DiY', DiY, 'DfY', DfY);

l.sw = load('f_fitted_off.mat');
l.sw = load('f_fitted_on.mat');
l.sw.Ronp = 112e-3;
l.sw.Rons = 112e-3;
% l.sw.Ronp = 90e-3;
% l.sw.Rons = 90e-3;

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



% l.tr.kc = 1.950522119690475;
l.tr.kc = 2.925783179535760;
l.tr.a = 1.448345520625370;
l.tr.b = 2.729356903304897;


% l.tr.a = 1.585;
% l.tr.b = 2.756;
% l.tr.kc = 0.603;
l.tr.int_ki = integral(@(theta) abs(cos(theta)).^l.tr.a,0,2*pi);
l.tr.ki = l.tr.kc/(2^(l.tr.b-l.tr.a)*(2*pi)^(l.tr.a-1)*l.tr.int_ki);
l.tr.Ac = 683e-6;
l.tr.Ve = 1;
% l.tr.Ve = 146.79e-6;
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


l.pr.dt = 300e-9;
l.pr.phi = deg2rad(67);
l.pr.Vi = 400;
l.pr.d = 300/400;
l.pr.fs = 100e3;
l.pr.Ldab = 61.5e-6;
l.pr.Ld1 = 1.4e-6;
l.pr.n = l.tr.Ns/l.tr.Np;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 700e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);

%% TRAFOOOOOOOOOOOOOO
trafoo = 'YY';

%%


color1 = [0.045 0.245 0.745]; % blue
color2 = [0.635 0.635 0.635]; % gray

if (trafoo == "YY")
    out = f_equationsYY(l);
    C = num2cell(out);
    [hbrm,HBrm,Ip,Is,iiRMS,iiME,ioRMS,ioME,Pm,idrm,ilrm,iLrm,iSwPrm,iSwSrm,P_core_tr,Bpk_tr,P_core_L,Bpk_L] = C{:};
    
    out = f_equationsYY_dt(l);
    
elseif (trafoo == "YD")
    out = f_equationsYD_dt(l);
elseif (trafoo == "DfD")
    out = f_equationsDfD_dt(l);
elseif (trafoo == "DiD")
    out = f_equationsDiD_dt(l);
elseif (trafoo == "DiY")
    out = f_equationsDiY_dt(l);
elseif (trafoo == "DfY")
    out = f_equationsDfY_dt(l);
end

[efc,ptot,Pm,cSw_p,cSw_s,sSw_p,sSw_s,Pil_cu,PiL_cu,PiLd_cu,Psc,P_core_L,P_core_tr] = f_loss(trafoo, out, l)


pse = l.tr.Ve*l.tr.kc*(l.pr.fs^l.tr.a)*((Bpk_tr)^l.tr.b)