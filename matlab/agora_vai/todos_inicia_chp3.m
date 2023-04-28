clear; close all; clc;

% as 6 opcoes estudas sao 
% 1: YY e iDY
% 2: YY e oDY
% 3: iDY e oDY 
% 4: YD e iDD
% 5: YD e oDD
% 6: iDD e oDD

syms Ld1 Ld2 n Lm Po t L1 L2 Ldab M fs Vi d dt real positive
syms phi real

color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue

pr = 200; %passo para plot

Io_limit = 10;
Vo_max = 500;
vec_vo_10 = 1:0.1:Vo_max;
vec_po_10 = Io_limit*vec_vo_10;
n_limit_10 = length(vec_vo_10);
n_10 = round(n_limit_10/4);


% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

addpath('utils')
addpath('utils_transf')


Vo = Vi*d;
Ts = 1/fs;

%% cria tudo do YY
%cria intervalo de angulo
YY.intervalo.f1 = [0 deg2rad(60)];
YY.intervalo.f2 = [deg2rad(60) deg2rad(89.999999)];

YY.trafo = 'YY'; %nome
%pega as funcoes analiticas
YY.phi_num = [sum(YY.intervalo.f1)/2,sum(YY.intervalo.f2)/2];  %[MUDAR]
[YY.x0s.f1, YY.ts.f1, YY.idab.f1, YY.hb.f1, YY.HB.f1, YY.Ip.f1, YY.Is.f1, YY.iME.f1, YY.idrm.f1, YY.ilrm.f1, YY.iLrm.f1, YY.iSwPrm.f1, YY.iSwSrm.f1] = simplify_YY(YY.phi_num(1));
[YY.x0s.f2, YY.ts.f2, YY.idab.f2, YY.hb.f2, YY.HB.f2, YY.Ip.f2, YY.Is.f2, YY.iME.f2, YY.idrm.f2, YY.ilrm.f2, YY.iLrm.f2, YY.iSwPrm.f2, YY.iSwSrm.f2] = simplify_YY(YY.phi_num(2));

%vetor do angulo
YY.vec_ph.f2 = min(YY.intervalo.f2):(max(YY.intervalo.f2)-min(YY.intervalo.f2))/pr:max(YY.intervalo.f2);
YY.vec_ph.f1 = min(YY.intervalo.f1):(max(YY.intervalo.f1)-min(YY.intervalo.f1))/pr:max(YY.intervalo.f1);
%corrente dos estados
YY.fx0s.f1 = matlabFunction(YY.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fx0s.f2 = matlabFunction(YY.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%tempo dos estados
YY.fts.f1 = matlabFunction(YY.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fts.f2 = matlabFunction(YY.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%equacoes da potencia no limite do zvs
YY.Ip_eq.f1 = YY.Ip.f1 == 0;
YY.Is_eq.f1 = YY.Is.f1 == 0;
YY.pot_eq.f1 = YY.iME.f1*Vo == Po;
YY.fpot_eq.f1 = matlabFunction(solve(YY.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIp_eq.f1 = matlabFunction(solve(YY.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIs_eq.f1 = matlabFunction(solve(YY.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.Ip_eq.f2 = YY.Ip.f2 == 0;
YY.Is_eq.f2 = YY.Is.f2 == 0;
YY.pot_eq.f2 = YY.iME.f2*Vo == Po;
YY.fpot_eq.f2 = matlabFunction(solve(YY.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIp_eq.f2 = matlabFunction(solve(YY.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIs_eq.f2 = matlabFunction(solve(YY.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%% cria tudo do iDY
%cria intervalo de angulo
iDY.intervalo.f1 = [deg2rad(-30) deg2rad(0)];
iDY.intervalo.f2 = [deg2rad(0) deg2rad(59.999999)];

iDY.trafo = 'iDY'; %nome
%pega as funcoes analiticas
iDY.phi_num = [sum(iDY.intervalo.f1)/2,sum(iDY.intervalo.f2)/2];  %[MUDAR]
[iDY.x0s.f1, iDY.ts.f1, iDY.idab.f1, iDY.hb.f1, iDY.HB.f1, iDY.Ip.f1, iDY.Is.f1, iDY.iME.f1, iDY.idrm.f1, iDY.ilrm.f1, iDY.iLrm.f1, iDY.iSwPrm.f1, iDY.iSwSrm.f1] = simplify_DiY(iDY.phi_num(1));
[iDY.x0s.f2, iDY.ts.f2, iDY.idab.f2, iDY.hb.f2, iDY.HB.f2, iDY.Ip.f2, iDY.Is.f2, iDY.iME.f2, iDY.idrm.f2, iDY.ilrm.f2, iDY.iLrm.f2, iDY.iSwPrm.f2, iDY.iSwSrm.f2] = simplify_DiY(iDY.phi_num(2));

%vetor do angulo
iDY.vec_ph.f2 = min(iDY.intervalo.f2):(max(iDY.intervalo.f2)-min(iDY.intervalo.f2))/pr:max(iDY.intervalo.f2);
iDY.vec_ph.f1 = min(iDY.intervalo.f1):(max(iDY.intervalo.f1)-min(iDY.intervalo.f1))/pr:max(iDY.intervalo.f1);
%corrente dos estados
iDY.fx0s.f1 = matlabFunction(iDY.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.fx0s.f2 = matlabFunction(iDY.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%tempo dos estados
iDY.fts.f1 = matlabFunction(iDY.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.fts.f2 = matlabFunction(iDY.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%equacoes da potencia no limite do zvs
iDY.Ip_eq.f1 = iDY.Ip.f1 == 0;
iDY.Is_eq.f1 = iDY.Is.f1 == 0;
iDY.pot_eq.f1 = iDY.iME.f1*Vo == Po;
iDY.fpot_eq.f1 = matlabFunction(solve(iDY.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.fIp_eq.f1 = matlabFunction(solve(iDY.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.fIs_eq.f1 = matlabFunction(solve(iDY.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.Ip_eq.f2 = iDY.Ip.f2 == 0;
iDY.Is_eq.f2 = iDY.Is.f2 == 0;
iDY.pot_eq.f2 = iDY.iME.f2*Vo == Po;
iDY.fpot_eq.f2 = matlabFunction(solve(iDY.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.fIp_eq.f2 = matlabFunction(solve(iDY.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.fIs_eq.f2 = matlabFunction(solve(iDY.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%% cria tudo do iDD
%cria intervalo de angulo
iDD.intervalo.f1 = [deg2rad(0) deg2rad(60)];
iDD.intervalo.f2 = [deg2rad(60) deg2rad(89.999999)];

iDD.trafo = 'iDD'; %nome
%pega as funcoes analiticas
iDD.phi_num = [sum(iDD.intervalo.f1)/2,sum(iDD.intervalo.f2)/2];  %[MUDAR]
[iDD.x0s.f1, iDD.ts.f1, iDD.idab.f1, iDD.hb.f1, iDD.HB.f1, iDD.Ip.f1, iDD.Is.f1, iDD.iME.f1, iDD.idrm.f1, iDD.ilrm.f1, iDD.iLrm.f1, iDD.iSwPrm.f1, iDD.iSwSrm.f1] = simplify_DiD(iDD.phi_num(1));
[iDD.x0s.f2, iDD.ts.f2, iDD.idab.f2, iDD.hb.f2, iDD.HB.f2, iDD.Ip.f2, iDD.Is.f2, iDD.iME.f2, iDD.idrm.f2, iDD.ilrm.f2, iDD.iLrm.f2, iDD.iSwPrm.f2, iDD.iSwSrm.f2] = simplify_DiD(iDD.phi_num(2));

%vetor do angulo
iDD.vec_ph.f2 = min(iDD.intervalo.f2):(max(iDD.intervalo.f2)-min(iDD.intervalo.f2))/pr:max(iDD.intervalo.f2);
iDD.vec_ph.f1 = min(iDD.intervalo.f1):(max(iDD.intervalo.f1)-min(iDD.intervalo.f1))/pr:max(iDD.intervalo.f1);
%corrente dos estados
iDD.fx0s.f1 = matlabFunction(iDD.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDD.fx0s.f2 = matlabFunction(iDD.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%tempo dos estados
iDD.fts.f1 = matlabFunction(iDD.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDD.fts.f2 = matlabFunction(iDD.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%equacoes da potencia no limite do zvs
iDD.Ip_eq.f1 = iDD.Ip.f1 == 0;
iDD.Is_eq.f1 = iDD.Is.f1 == 0;
iDD.pot_eq.f1 = iDD.iME.f1*Vo == Po;
iDD.fpot_eq.f1 = matlabFunction(solve(iDD.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDD.fIp_eq.f1 = matlabFunction(solve(iDD.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDD.fIs_eq.f1 = matlabFunction(solve(iDD.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDD.Ip_eq.f2 = iDD.Ip.f2 == 0;
iDD.Is_eq.f2 = iDD.Is.f2 == 0;
iDD.pot_eq.f2 = iDD.iME.f2*Vo == Po;
iDD.fpot_eq.f2 = matlabFunction(solve(iDD.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDD.fIp_eq.f2 = matlabFunction(solve(iDD.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDD.fIs_eq.f2 = matlabFunction(solve(iDD.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%% cria tudo do oDD
%cria intervalo de angulo
oDD.intervalo.f1 = [deg2rad(0) deg2rad(60)];
oDD.intervalo.f2 = [deg2rad(60) deg2rad(89.999999)];

oDD.trafo = 'oDD'; %nome
%pega as funcoes analiticas
oDD.phi_num = [sum(oDD.intervalo.f1)/2,sum(oDD.intervalo.f2)/2];  %[MUDAR]
[oDD.x0s.f1, oDD.ts.f1, oDD.idab.f1, oDD.hb.f1, oDD.HB.f1, oDD.Ip.f1, oDD.Is.f1, oDD.iME.f1, oDD.idrm.f1, oDD.ilrm.f1, oDD.iLrm.f1, oDD.iSwPrm.f1, oDD.iSwSrm.f1] = simplify_DfD(oDD.phi_num(1));
[oDD.x0s.f2, oDD.ts.f2, oDD.idab.f2, oDD.hb.f2, oDD.HB.f2, oDD.Ip.f2, oDD.Is.f2, oDD.iME.f2, oDD.idrm.f2, oDD.ilrm.f2, oDD.iLrm.f2, oDD.iSwPrm.f2, oDD.iSwSrm.f2] = simplify_DfD(oDD.phi_num(2));

%vetor do angulo
oDD.vec_ph.f2 = min(oDD.intervalo.f2):(max(oDD.intervalo.f2)-min(oDD.intervalo.f2))/pr:max(oDD.intervalo.f2);
oDD.vec_ph.f1 = min(oDD.intervalo.f1):(max(oDD.intervalo.f1)-min(oDD.intervalo.f1))/pr:max(oDD.intervalo.f1);
%corrente dos estados
oDD.fx0s.f1 = matlabFunction(oDD.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDD.fx0s.f2 = matlabFunction(oDD.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%tempo dos estados
oDD.fts.f1 = matlabFunction(oDD.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDD.fts.f2 = matlabFunction(oDD.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%equacoes da potencia no limite do zvs
oDD.Ip_eq.f1 = oDD.Ip.f1 == 0;
oDD.Is_eq.f1 = oDD.Is.f1 == 0;
oDD.pot_eq.f1 = oDD.iME.f1*Vo == Po;
oDD.fpot_eq.f1 = matlabFunction(solve(oDD.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDD.fIp_eq.f1 = matlabFunction(solve(oDD.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDD.fIs_eq.f1 = matlabFunction(solve(oDD.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDD.Ip_eq.f2 = oDD.Ip.f2 == 0;
oDD.Is_eq.f2 = oDD.Is.f2 == 0;
oDD.pot_eq.f2 = oDD.iME.f2*Vo == Po;
oDD.fpot_eq.f2 = matlabFunction(solve(oDD.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDD.fIp_eq.f2 = matlabFunction(solve(oDD.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDD.fIs_eq.f2 = matlabFunction(solve(oDD.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%% cria tudo do YD
%cria intervalo de angulo
YD.intervalo.f1 = [deg2rad(30) deg2rad(60)];
YD.intervalo.f2 = [deg2rad(60) deg2rad(119.999999)];

YD.trafo = 'YD'; %nome
%pega as funcoes analiticas
YD.phi_num = [sum(YD.intervalo.f1)/2,sum(YD.intervalo.f2)/2];  %[MUDAR]
[YD.x0s.f1, YD.ts.f1, YD.idab.f1, YD.hb.f1, YD.HB.f1, YD.Ip.f1, YD.Is.f1, YD.iME.f1, YD.idrm.f1, YD.ilrm.f1, YD.iLrm.f1, YD.iSwPrm.f1, YD.iSwSrm.f1] = simplify_YD(YD.phi_num(1));
[YD.x0s.f2, YD.ts.f2, YD.idab.f2, YD.hb.f2, YD.HB.f2, YD.Ip.f2, YD.Is.f2, YD.iME.f2, YD.idrm.f2, YD.ilrm.f2, YD.iLrm.f2, YD.iSwPrm.f2, YD.iSwSrm.f2] = simplify_YD(YD.phi_num(2));

%vetor do angulo
YD.vec_ph.f2 = min(YD.intervalo.f2):(max(YD.intervalo.f2)-min(YD.intervalo.f2))/pr:max(YD.intervalo.f2);
YD.vec_ph.f1 = min(YD.intervalo.f1):(max(YD.intervalo.f1)-min(YD.intervalo.f1))/pr:max(YD.intervalo.f1);
%corrente dos estados
YD.fx0s.f1 = matlabFunction(YD.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fx0s.f2 = matlabFunction(YD.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%tempo dos estados
YD.fts.f1 = matlabFunction(YD.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fts.f2 = matlabFunction(YD.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%equacoes da potencia no limite do zvs
YD.Ip_eq.f1 = YD.Ip.f1 == 0;
YD.Is_eq.f1 = YD.Is.f1 == 0;
YD.pot_eq.f1 = YD.iME.f1*Vo == Po;
YD.fpot_eq.f1 = matlabFunction(solve(YD.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIp_eq.f1 = matlabFunction(solve(YD.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIs_eq.f1 = matlabFunction(solve(YD.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.Ip_eq.f2 = YD.Ip.f2 == 0;
YD.Is_eq.f2 = YD.Is.f2 == 0;
YD.pot_eq.f2 = YD.iME.f2*Vo == Po;
YD.fpot_eq.f2 = matlabFunction(solve(YD.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIp_eq.f2 = matlabFunction(solve(YD.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIs_eq.f2 = matlabFunction(solve(YD.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%% cria tudo do oDY
%cria intervalo de angulo
oDY.intervalo.f1 = [deg2rad(-30) deg2rad(0)];
oDY.intervalo.f2 = [deg2rad(0) deg2rad(59.999999)];

oDY.trafo = 'oDY'; %nome
%pega as funcoes analiticas
oDY.phi_num = [sum(oDY.intervalo.f1)/2,sum(oDY.intervalo.f2)/2];  %[MUDAR]
[oDY.x0s.f1, oDY.ts.f1, oDY.idab.f1, oDY.hb.f1, oDY.HB.f1, oDY.Ip.f1, oDY.Is.f1, oDY.iME.f1, oDY.idrm.f1, oDY.ilrm.f1, oDY.iLrm.f1, oDY.iSwPrm.f1, oDY.iSwSrm.f1] = simplify_DfY(oDY.phi_num(1));
[oDY.x0s.f2, oDY.ts.f2, oDY.idab.f2, oDY.hb.f2, oDY.HB.f2, oDY.Ip.f2, oDY.Is.f2, oDY.iME.f2, oDY.idrm.f2, oDY.ilrm.f2, oDY.iLrm.f2, oDY.iSwPrm.f2, oDY.iSwSrm.f2] = simplify_DfY(oDY.phi_num(2));

%vetor do angulo
oDY.vec_ph.f2 = min(oDY.intervalo.f2):(max(oDY.intervalo.f2)-min(oDY.intervalo.f2))/pr:max(oDY.intervalo.f2);
oDY.vec_ph.f1 = min(oDY.intervalo.f1):(max(oDY.intervalo.f1)-min(oDY.intervalo.f1))/pr:max(oDY.intervalo.f1);
%corrente dos estados
oDY.fx0s.f1 = matlabFunction(oDY.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDY.fx0s.f2 = matlabFunction(oDY.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%tempo dos estados
oDY.fts.f1 = matlabFunction(oDY.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDY.fts.f2 = matlabFunction(oDY.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
%equacoes da potencia no limite do zvs
oDY.Ip_eq.f1 = oDY.Ip.f1 == 0;
oDY.Is_eq.f1 = oDY.Is.f1 == 0;
oDY.pot_eq.f1 = oDY.iME.f1*Vo == Po;
oDY.fpot_eq.f1 = matlabFunction(solve(oDY.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDY.fIp_eq.f1 = matlabFunction(solve(oDY.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDY.fIs_eq.f1 = matlabFunction(solve(oDY.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDY.Ip_eq.f2 = oDY.Ip.f2 == 0;
oDY.Is_eq.f2 = oDY.Is.f2 == 0;
oDY.pot_eq.f2 = oDY.iME.f2*Vo == Po;
oDY.fpot_eq.f2 = matlabFunction(solve(oDY.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDY.fIp_eq.f2 = matlabFunction(solve(oDY.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
oDY.fIs_eq.f2 = matlabFunction(solve(oDY.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});



%% opcao 1: YY e iDY

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 60e-6;
Ld1_num = 1.4e-6;
n_num = 5/9;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/5);

figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
fill([YY.lim_pot.f1 fliplr(YY.lim_pot.f1)], [Vi_num*YY.lim_p.f1 fliplr(Vi_num*YY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
fill([YY.lim_pot.f2 fliplr(YY.lim_pot.f2)], [Vi_num*YY.lim_p.f2 fliplr(Vi_num*YY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(YY.lim_pot.f1,Vi_num*YY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f1,Vi_num*YY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(YY.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)
    
iDY.lim_p.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_s.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_pot.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
fill([iDY.lim_pot.f1 fliplr(iDY.lim_pot.f1)], [Vi_num*iDY.lim_p.f1 fliplr(Vi_num*iDY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
iDY.lim_p.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_s.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_pot.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
fill([iDY.lim_pot.f2 fliplr(iDY.lim_pot.f2)], [Vi_num*iDY.lim_p.f2 fliplr(Vi_num*iDY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(iDY.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'YY','iDY'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',YY.trafo,'_',iDY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% opcao 2: YY e oDY

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 60e-6;
Ld1_num = 1.4e-6;
n_num = 5/9;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/6);

figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
fill([YY.lim_pot.f1 fliplr(YY.lim_pot.f1)], [Vi_num*YY.lim_p.f1 fliplr(Vi_num*YY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
fill([YY.lim_pot.f2 fliplr(YY.lim_pot.f2)], [Vi_num*YY.lim_p.f2 fliplr(Vi_num*YY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(YY.lim_pot.f1,Vi_num*YY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f1,Vi_num*YY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(YY.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)

oDY.lim_p.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
oDY.lim_s.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
oDY.lim_pot.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
fill([oDY.lim_pot.f1 fliplr(oDY.lim_pot.f1)], [Vi_num*oDY.lim_p.f1 fliplr(Vi_num*oDY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
oDY.lim_p.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
oDY.lim_s.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
oDY.lim_pot.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
fill([oDY.lim_pot.f2 fliplr(oDY.lim_pot.f2)], [Vi_num*oDY.lim_p.f2 fliplr(Vi_num*oDY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(oDY.lim_pot.f1,Vi_num*oDY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f2,Vi_num*oDY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f1,Vi_num*oDY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f2,Vi_num*oDY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(oDY.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'YY','oDY'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',YY.trafo,'_',oDY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% opcao 3: iDY e oDY

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 90e-6;
Ld1_num = 1.4e-6;
n_num = 5/9;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/2.5);

figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

iDY.lim_p.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_s.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_pot.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
fill([iDY.lim_pot.f1 fliplr(iDY.lim_pot.f1)], [Vi_num*iDY.lim_p.f1 fliplr(Vi_num*iDY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
iDY.lim_p.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_s.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_pot.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
fill([iDY.lim_pot.f2 fliplr(iDY.lim_pot.f2)], [Vi_num*iDY.lim_p.f2 fliplr(Vi_num*iDY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(iDY.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)

oDY.lim_p.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
oDY.lim_s.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
oDY.lim_pot.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
fill([oDY.lim_pot.f1 fliplr(oDY.lim_pot.f1)], [Vi_num*oDY.lim_p.f1 fliplr(Vi_num*oDY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
oDY.lim_p.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
oDY.lim_s.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
oDY.lim_pot.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
fill([oDY.lim_pot.f2 fliplr(oDY.lim_pot.f2)], [Vi_num*oDY.lim_p.f2 fliplr(Vi_num*oDY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(oDY.lim_pot.f1,Vi_num*oDY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f2,Vi_num*oDY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f1,Vi_num*oDY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f2,Vi_num*oDY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(oDY.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'iDY','oDY'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',iDY.trafo,'_',oDY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% opcao 4: YD e iDD

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 100e-6;
Ld1_num = 1.4e-6;
n_num = 1.05;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/6);

figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

YD.lim_p.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
YD.lim_s.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
YD.lim_pot.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
fill([YD.lim_pot.f1 fliplr(YD.lim_pot.f1)], [Vi_num*YD.lim_p.f1 fliplr(Vi_num*YD.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
YD.lim_p.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
YD.lim_s.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
YD.lim_pot.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
fill([YD.lim_pot.f2 fliplr(YD.lim_pot.f2)], [Vi_num*YD.lim_p.f2 fliplr(Vi_num*YD.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(YD.lim_pot.f1,Vi_num*YD.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YD.lim_pot.f2,Vi_num*YD.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YD.lim_pot.f1,Vi_num*YD.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YD.lim_pot.f2,Vi_num*YD.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(YD.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)

iDD.lim_p.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
iDD.lim_s.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
iDD.lim_pot.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
fill([iDD.lim_pot.f1 fliplr(iDD.lim_pot.f1)], [Vi_num*iDD.lim_p.f1 fliplr(Vi_num*iDD.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
iDD.lim_p.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
iDD.lim_s.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
iDD.lim_pot.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
fill([iDD.lim_pot.f2 fliplr(iDD.lim_pot.f2)], [Vi_num*iDD.lim_p.f2 fliplr(Vi_num*iDD.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(iDD.lim_pot.f1,Vi_num*iDD.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDD.lim_pot.f2,Vi_num*iDD.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDD.lim_pot.f1,Vi_num*iDD.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDD.lim_pot.f2,Vi_num*iDD.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(iDD.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'YD','iDD'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',YD.trafo,'_',iDD.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% opcao 5: YD e oDD

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 60e-6;
Ld1_num = 1.4e-6;
n_num = 1;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/6);

figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

YD.lim_p.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
YD.lim_s.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
YD.lim_pot.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
fill([YD.lim_pot.f1 fliplr(YD.lim_pot.f1)], [Vi_num*YD.lim_p.f1 fliplr(Vi_num*YD.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
YD.lim_p.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
YD.lim_s.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
YD.lim_pot.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
fill([YD.lim_pot.f2 fliplr(YD.lim_pot.f2)], [Vi_num*YD.lim_p.f2 fliplr(Vi_num*YD.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(YD.lim_pot.f1,Vi_num*YD.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YD.lim_pot.f2,Vi_num*YD.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YD.lim_pot.f1,Vi_num*YD.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YD.lim_pot.f2,Vi_num*YD.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(YD.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)

oDD.lim_p.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
oDD.lim_s.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
oDD.lim_pot.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
fill([oDD.lim_pot.f1 fliplr(oDD.lim_pot.f1)], [Vi_num*oDD.lim_p.f1 fliplr(Vi_num*oDD.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
oDD.lim_p.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
oDD.lim_s.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
oDD.lim_pot.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
fill([oDD.lim_pot.f2 fliplr(oDD.lim_pot.f2)], [Vi_num*oDD.lim_p.f2 fliplr(Vi_num*oDD.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(oDD.lim_pot.f1,Vi_num*oDD.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f2,Vi_num*oDD.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f1,Vi_num*oDD.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f2,Vi_num*oDD.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(oDD.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'YD','oDD'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',YD.trafo,'_',oDD.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% opcao 6: iDD e oDD

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 100e-6;
Ld1_num = 1.4e-6;
n_num = 1;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/2.8);

figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

iDD.lim_p.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
iDD.lim_s.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
iDD.lim_pot.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
fill([iDD.lim_pot.f1 fliplr(iDD.lim_pot.f1)], [Vi_num*iDD.lim_p.f1 fliplr(Vi_num*iDD.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
iDD.lim_p.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
iDD.lim_s.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
iDD.lim_pot.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
fill([iDD.lim_pot.f2 fliplr(iDD.lim_pot.f2)], [Vi_num*iDD.lim_p.f2 fliplr(Vi_num*iDD.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(iDD.lim_pot.f1,Vi_num*iDD.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDD.lim_pot.f2,Vi_num*iDD.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDD.lim_pot.f1,Vi_num*iDD.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDD.lim_pot.f2,Vi_num*iDD.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(iDD.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)

oDD.lim_p.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
oDD.lim_s.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
oDD.lim_pot.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
fill([oDD.lim_pot.f1 fliplr(oDD.lim_pot.f1)], [Vi_num*oDD.lim_p.f1 fliplr(Vi_num*oDD.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
oDD.lim_p.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
oDD.lim_s.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
oDD.lim_pot.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
fill([oDD.lim_pot.f2 fliplr(oDD.lim_pot.f2)], [Vi_num*oDD.lim_p.f2 fliplr(Vi_num*oDD.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(oDD.lim_pot.f1,Vi_num*oDD.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f2,Vi_num*oDD.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f1,Vi_num*oDD.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f2,Vi_num*oDD.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(oDD.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'iDD','oDD'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',iDD.trafo,'_',oDD.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% opcao 7: iDD e oDD e YD

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 90e-6;
Ld1_num = 1.4e-6;
n_num = 1.05;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/10);

figure
hold on

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:3
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

iDD.lim_p.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
iDD.lim_s.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
iDD.lim_pot.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f1);
fill([iDD.lim_pot.f1 fliplr(iDD.lim_pot.f1)], [Vi_num*iDD.lim_p.f1 fliplr(Vi_num*iDD.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
iDD.lim_p.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
iDD.lim_s.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
iDD.lim_pot.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDD.vec_ph.f2);
fill([iDD.lim_pot.f2 fliplr(iDD.lim_pot.f2)], [Vi_num*iDD.lim_p.f2 fliplr(Vi_num*iDD.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(iDD.lim_pot.f1,Vi_num*iDD.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDD.lim_pot.f2,Vi_num*iDD.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDD.lim_pot.f1,Vi_num*iDD.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDD.lim_pot.f2,Vi_num*iDD.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(iDD.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)

oDD.lim_p.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
oDD.lim_s.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
oDD.lim_pot.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f1);
fill([oDD.lim_pot.f1 fliplr(oDD.lim_pot.f1)], [Vi_num*oDD.lim_p.f1 fliplr(Vi_num*oDD.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
oDD.lim_p.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
oDD.lim_s.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
oDD.lim_pot.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDD.vec_ph.f2);
fill([oDD.lim_pot.f2 fliplr(oDD.lim_pot.f2)], [Vi_num*oDD.lim_p.f2 fliplr(Vi_num*oDD.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(oDD.lim_pot.f1,Vi_num*oDD.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f2,Vi_num*oDD.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f1,Vi_num*oDD.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDD.lim_pot.f2,Vi_num*oDD.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(oDD.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

YD.lim_p.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
YD.lim_s.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
YD.lim_pot.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f1);
fill([YD.lim_pot.f1 fliplr(YD.lim_pot.f1)], [Vi_num*YD.lim_p.f1 fliplr(Vi_num*YD.lim_s.f1)],jetcustom(3,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
YD.lim_p.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
YD.lim_s.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
YD.lim_pot.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YD.vec_ph.f2);
fill([YD.lim_pot.f2 fliplr(YD.lim_pot.f2)], [Vi_num*YD.lim_p.f2 fliplr(Vi_num*YD.lim_s.f2)],jetcustom(3,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(YD.lim_pot.f1,Vi_num*YD.lim_p.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(YD.lim_pot.f2,Vi_num*YD.lim_p.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(YD.lim_pot.f1,Vi_num*YD.lim_s.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(YD.lim_pot.f2,Vi_num*YD.lim_s.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
xline(max(YD.lim_pot.f2),'Color',jetcustom(3,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'iDD','oDD','YD'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',iDD.trafo,'_',oDD.trafo,'_',YD.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% opcao 8: iDY e oDY e YY

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 90e-6;
Ld1_num = 1.4e-6;
n_num = 5/9;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

n_10 = round(n_limit_10/10);

figure
hold on

for ii=1:3
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

iDY.lim_p.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_s.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_pot.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
fill([iDY.lim_pot.f1 fliplr(iDY.lim_pot.f1)], [Vi_num*iDY.lim_p.f1 fliplr(Vi_num*iDY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
iDY.lim_p.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_s.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_pot.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
fill([iDY.lim_pot.f2 fliplr(iDY.lim_pot.f2)], [Vi_num*iDY.lim_p.f2 fliplr(Vi_num*iDY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(iDY.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)

oDY.lim_p.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
oDY.lim_s.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
oDY.lim_pot.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f1);
fill([oDY.lim_pot.f1 fliplr(oDY.lim_pot.f1)], [Vi_num*oDY.lim_p.f1 fliplr(Vi_num*oDY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
oDY.lim_p.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
oDY.lim_s.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
oDY.lim_pot.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,oDY.vec_ph.f2);
fill([oDY.lim_pot.f2 fliplr(oDY.lim_pot.f2)], [Vi_num*oDY.lim_p.f2 fliplr(Vi_num*oDY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(oDY.lim_pot.f1,Vi_num*oDY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f2,Vi_num*oDY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f1,Vi_num*oDY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(oDY.lim_pot.f2,Vi_num*oDY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(oDY.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
fill([YY.lim_pot.f1 fliplr(YY.lim_pot.f1)], [Vi_num*YY.lim_p.f1 fliplr(Vi_num*YY.lim_s.f1)],jetcustom(3,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
fill([YY.lim_pot.f2 fliplr(YY.lim_pot.f2)], [Vi_num*YY.lim_p.f2 fliplr(Vi_num*YY.lim_s.f2)],jetcustom(3,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(YY.lim_pot.f1,Vi_num*YY.lim_p.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_p.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(YY.lim_pot.f1,Vi_num*YY.lim_s.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_s.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
xline(max(YY.lim_pot.f2),'Color',jetcustom(3,:),'LineWidth',1.5)

plot(vec_po_10, vec_vo_10,'LineWidth',1.5,'Color',[0 0 0])
text(vec_po_10(n_10)+100,vec_vo_10(n_10),'$\leftarrow\,10\,$[A]','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 4500])
legend({'iDD','oDD','YD'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\',iDY.trafo,'_',oDY.trafo,'_',YY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');


%% analisar eles sozinhos, parametros todos iguais

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 5e-6;
n_num = 1;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = [400e-6, 1e-3 10e-3];
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

%% YY sozinho
figure
hold on
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
    YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
    YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
    fill([YY.lim_pot.f1 fliplr(YY.lim_pot.f1)], [Vi_num*YY.lim_p.f1 fliplr(Vi_num*YY.lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
    YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
    YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
    fill([YY.lim_pot.f2 fliplr(YY.lim_pot.f2)], [Vi_num*YY.lim_p.f2 fliplr(Vi_num*YY.lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(YY.lim_pot.f1,Vi_num*YY.lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(YY.lim_pot.f2,Vi_num*YY.lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(YY.lim_pot.f1,Vi_num*YY.lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(YY.lim_pot.f2,Vi_num*YY.lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end
hold off
ylim([0 1.5*Vi_num])
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

text(200,550,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(200,220,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(1000,400,'SS in both','Interpreter', 'Latex','FontSize', 14) 

file_name = append('figure\comparacoes\sozinho_',YY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');


%% iDY sozinho
figure
hold on
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    iDY.lim_p.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDY.vec_ph.f1);
    iDY.lim_s.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDY.vec_ph.f1);
    iDY.lim_pot.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDY.vec_ph.f1);
    fill([iDY.lim_pot.f1 fliplr(iDY.lim_pot.f1)], [Vi_num*iDY.lim_p.f1 fliplr(Vi_num*iDY.lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    iDY.lim_p.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDY.vec_ph.f2);
    iDY.lim_s.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDY.vec_ph.f2);
    iDY.lim_pot.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDY.vec_ph.f2);
    fill([iDY.lim_pot.f2 fliplr(iDY.lim_pot.f2)], [Vi_num*iDY.lim_p.f2 fliplr(Vi_num*iDY.lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(iDY.lim_pot.f1,Vi_num*iDY.lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(iDY.lim_pot.f2,Vi_num*iDY.lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(iDY.lim_pot.f1,Vi_num*iDY.lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(iDY.lim_pot.f2,Vi_num*iDY.lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end
hold off
ylim([0 1000])
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

text(300,900,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(300,400,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(2500,700,'SS in both','Interpreter', 'Latex','FontSize', 14) 

file_name = append('figure\comparacoes\sozinho_',iDY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% oDY sozinho
figure
hold on
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    oDY.lim_p.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDY.vec_ph.f1);
    oDY.lim_s.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDY.vec_ph.f1);
    oDY.lim_pot.f1 = ones(1,length(oDY.vec_ph.f1)).*oDY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDY.vec_ph.f1);
    fill([oDY.lim_pot.f1 fliplr(oDY.lim_pot.f1)], [Vi_num*oDY.lim_p.f1 fliplr(Vi_num*oDY.lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    oDY.lim_p.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDY.vec_ph.f2);
    oDY.lim_s.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDY.vec_ph.f2);
    oDY.lim_pot.f2 = ones(1,length(oDY.vec_ph.f2)).*oDY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDY.vec_ph.f2);
    fill([oDY.lim_pot.f2 fliplr(oDY.lim_pot.f2)], [Vi_num*oDY.lim_p.f2 fliplr(Vi_num*oDY.lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(oDY.lim_pot.f1,Vi_num*oDY.lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(oDY.lim_pot.f2,Vi_num*oDY.lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(oDY.lim_pot.f1,Vi_num*oDY.lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(oDY.lim_pot.f2,Vi_num*oDY.lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end
hold off
ylim([0 1000])
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

text(200,900,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(200,350,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(800,700,'SS in both','Interpreter', 'Latex','FontSize', 14) 

file_name = append('figure\comparacoes\sozinho_',oDY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% YD sozinho
figure
hold on
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    YD.lim_p.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YD.vec_ph.f1);
    YD.lim_s.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YD.vec_ph.f1);
    YD.lim_pot.f1 = ones(1,length(YD.vec_ph.f1)).*YD.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YD.vec_ph.f1);
    fill([YD.lim_pot.f1 fliplr(YD.lim_pot.f1)], [Vi_num*YD.lim_p.f1 fliplr(Vi_num*YD.lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    YD.lim_p.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YD.vec_ph.f2);
    YD.lim_s.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YD.vec_ph.f2);
    YD.lim_pot.f2 = ones(1,length(YD.vec_ph.f2)).*YD.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YD.vec_ph.f2);
    fill([YD.lim_pot.f2 fliplr(YD.lim_pot.f2)], [Vi_num*YD.lim_p.f2 fliplr(Vi_num*YD.lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(YD.lim_pot.f1,Vi_num*YD.lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(YD.lim_pot.f2,Vi_num*YD.lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(YD.lim_pot.f1,Vi_num*YD.lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(YD.lim_pot.f2,Vi_num*YD.lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end
hold off
ylim([0 400])
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

text(300,340,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(300,138,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(2500,250,'SS in both','Interpreter', 'Latex','FontSize', 14) 

file_name = append('figure\comparacoes\sozinho_',YD.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% iDD sozinho
figure
hold on
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    iDD.lim_p.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDD.vec_ph.f1);
    iDD.lim_s.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDD.vec_ph.f1);
    iDD.lim_pot.f1 = ones(1,length(iDD.vec_ph.f1)).*iDD.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDD.vec_ph.f1);
    fill([iDD.lim_pot.f1 fliplr(iDD.lim_pot.f1)], [Vi_num*iDD.lim_p.f1 fliplr(Vi_num*iDD.lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    iDD.lim_p.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDD.vec_ph.f2);
    iDD.lim_s.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDD.vec_ph.f2);
    iDD.lim_pot.f2 = ones(1,length(iDD.vec_ph.f2)).*iDD.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,iDD.vec_ph.f2);
    fill([iDD.lim_pot.f2 fliplr(iDD.lim_pot.f2)], [Vi_num*iDD.lim_p.f2 fliplr(Vi_num*iDD.lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(iDD.lim_pot.f1,Vi_num*iDD.lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(iDD.lim_pot.f2,Vi_num*iDD.lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(iDD.lim_pot.f1,Vi_num*iDD.lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(iDD.lim_pot.f2,Vi_num*iDD.lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end
hold off
ylim([0 1000])
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

text(800,600,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(800,150,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(4000,400,'SS in both','Interpreter', 'Latex','FontSize', 14) 

file_name = append('figure\comparacoes\sozinho_',iDD.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% oDD sozinho
figure
hold on
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    oDD.lim_p.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDD.vec_ph.f1);
    oDD.lim_s.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDD.vec_ph.f1);
    oDD.lim_pot.f1 = ones(1,length(oDD.vec_ph.f1)).*oDD.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDD.vec_ph.f1);
    fill([oDD.lim_pot.f1 fliplr(oDD.lim_pot.f1)], [Vi_num*oDD.lim_p.f1 fliplr(Vi_num*oDD.lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    oDD.lim_p.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDD.vec_ph.f2);
    oDD.lim_s.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDD.vec_ph.f2);
    oDD.lim_pot.f2 = ones(1,length(oDD.vec_ph.f2)).*oDD.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,oDD.vec_ph.f2);
    fill([oDD.lim_pot.f2 fliplr(oDD.lim_pot.f2)], [Vi_num*oDD.lim_p.f2 fliplr(Vi_num*oDD.lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(oDD.lim_pot.f1,Vi_num*oDD.lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(oDD.lim_pot.f2,Vi_num*oDD.lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(oDD.lim_pot.f1,Vi_num*oDD.lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(oDD.lim_pot.f2,Vi_num*oDD.lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end
hold off
ylim([0 1000])
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

text(300,600,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(300,120,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(1500,400,'SS in both','Interpreter', 'Latex','FontSize', 14) 

file_name = append('figure\comparacoes\sozinho_',oDD.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%%

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61.5e-6;
Ld1_num = 1.4e-6;
n_num = 5/9;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = .7e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);


%pts experimentais testados
pt_vo = [400, 300];
pt_io = [10,    10];
pt_po = pt_vo.*pt_io;
io_offset = 0.5;
po_offset = 200;

%% calcular phi teorico
syms Ptarget

iDY.phi_target.f1 = matlabFunction(solve(d*Vi*iDY.iME.f1 == Ptarget,phi), 'vars', {L1,L2,Ldab,M,Vi,d,fs,Ptarget});
iDY.phi_target.f2 = matlabFunction(solve(d*Vi*iDY.iME.f2 == Ptarget,phi), 'vars', {L1,L2,Ldab,M,Vi,d,fs,Ptarget});
YY.phi_target.f1 = matlabFunction(solve(d*Vi*YY.iME.f1 == Ptarget,phi), 'vars', {L1,L2,Ldab,M,Vi,d,fs,Ptarget});
YY.phi_target.f2 = matlabFunction(solve(d*Vi*YY.iME.f2 == Ptarget,phi), 'vars', {L1,L2,Ldab,M,Vi,d,fs,Ptarget});

target = 2;
opcoes = YY.phi_target.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,pt_vo(target)/400,fs_num,pt_po(target))*180/pi;
opcoes(1)
opcoes(2)
opcoes = YY.phi_target.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,pt_vo(target)/400,fs_num,pt_po(target))*180/pi;
opcoes(1)
opcoes(2)
%%

n_10 = round(n_limit_10/5);

figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
fill([YY.lim_pot.f1 fliplr(YY.lim_pot.f1)], [Vi_num*YY.lim_p.f1 fliplr(Vi_num*YY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
fill([YY.lim_pot.f2 fliplr(YY.lim_pot.f2)], [Vi_num*YY.lim_p.f2 fliplr(Vi_num*YY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
plot(YY.lim_pot.f1,Vi_num*YY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f1,Vi_num*YY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2,Vi_num*YY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
xline(max(YY.lim_pot.f2),'Color',jetcustom(1,:),'LineWidth',1.5)
    
iDY.lim_p.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_s.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_pot.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
fill([iDY.lim_pot.f1 fliplr(iDY.lim_pot.f1)], [Vi_num*iDY.lim_p.f1 fliplr(Vi_num*iDY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
iDY.lim_p.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_s.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_pot.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
fill([iDY.lim_pot.f2 fliplr(iDY.lim_pot.f2)], [Vi_num*iDY.lim_p.f2 fliplr(Vi_num*iDY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f1,Vi_num*iDY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2,Vi_num*iDY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
xline(max(iDY.lim_pot.f2),'Color',jetcustom(2,:),'LineWidth',1.5)

scatter(pt_po,pt_vo,'filled','square','MarkerEdgeColor','k','MarkerFaceColor','k')
text(pt_po(1)+po_offset,pt_vo(1),'$\leftarrow$ A','Interpreter', 'Latex','FontSize', 16) 
text(pt_po(2)+po_offset,pt_vo(2),'$\leftarrow$ B','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 5000])
legend({'YY','iDY'},'Location','southwest','FontSize', 14)
grid on
grid minor
xlabel('$P_o\,$[W]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\pontos_testados_Po_',YY.trafo,'_',iDY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');



%%
figure
hold on

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

for ii=1:2
    h = fill([1 1], [1 1],jetcustom(ii,:),'Edgecolor', 'none');
    h.FaceAlpha = 0.5;
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f1);
YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);
YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,YY.vec_ph.f2);

fill([YY.lim_pot.f1./(Vi_num*YY.lim_p.f1) fliplr(YY.lim_pot.f1./(Vi_num*YY.lim_s.f1))], [Vi_num*YY.lim_p.f1 fliplr(Vi_num*YY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
fill([YY.lim_pot.f2./(Vi_num*YY.lim_p.f2) fliplr(YY.lim_pot.f2./(Vi_num*YY.lim_s.f2))], [Vi_num*YY.lim_p.f2 fliplr(Vi_num*YY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(YY.lim_pot.f1./(Vi_num*YY.lim_p.f1),Vi_num*YY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2./(Vi_num*YY.lim_p.f2),Vi_num*YY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f1./(Vi_num*YY.lim_s.f1),Vi_num*YY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.lim_pot.f2./(Vi_num*YY.lim_s.f2),Vi_num*YY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
    
iDY.lim_p.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIp_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_s.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fIs_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_pot.f1 = ones(1,length(iDY.vec_ph.f1)).*iDY.fpot_eq.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f1);
iDY.lim_p.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIp_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_s.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fIs_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);
iDY.lim_pot.f2 = ones(1,length(iDY.vec_ph.f2)).*iDY.fpot_eq.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,iDY.vec_ph.f2);

fill([iDY.lim_pot.f1./(Vi_num*iDY.lim_p.f1) fliplr(iDY.lim_pot.f1./(Vi_num*iDY.lim_s.f1))], [Vi_num*iDY.lim_p.f1 fliplr(Vi_num*iDY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
fill([iDY.lim_pot.f2./(Vi_num*iDY.lim_p.f2) fliplr(iDY.lim_pot.f2./(Vi_num*iDY.lim_s.f2))], [Vi_num*iDY.lim_p.f2 fliplr(Vi_num*iDY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(iDY.lim_pot.f1./(Vi_num*iDY.lim_p.f1),Vi_num*iDY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2./(Vi_num*iDY.lim_p.f2),Vi_num*iDY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f1./(Vi_num*iDY.lim_s.f1),Vi_num*iDY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(iDY.lim_pot.f2./(Vi_num*iDY.lim_s.f2),Vi_num*iDY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)

scatter(pt_io,pt_vo,'filled','square','MarkerEdgeColor','k','MarkerFaceColor','k')
text(pt_io(1)+io_offset,pt_vo(1),'$\leftarrow$ A','Interpreter', 'Latex','FontSize', 16) 
text(pt_io(2)+io_offset,pt_vo(2),'$\leftarrow$ B','Interpreter', 'Latex','FontSize', 16) 
% text(pt_io(3)+io_offset,pt_vo(3),'$\leftarrow$ C','Interpreter', 'Latex','FontSize', 16) 

hold off
ylim([0 1.2*Vi_num])
xlim([0 12])
legend({'YY','iDY'},'Location','southeast','FontSize', 14)
grid on
grid minor
xlabel('$I_o\,$[A]')
ylabel('$V_o\,$[V]')
set(gca, 'FontSize', 20)

file_name = append('figure\comparacoes\pontos_testados_Io_',YY.trafo,'_',iDY.trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% potencia de saida em funcao do angulo
pr2 = pr*100;
iDY.vec_ph2.f2 = min(iDY.intervalo.f2):(max(iDY.intervalo.f2)-min(iDY.intervalo.f2))/pr2:max(iDY.intervalo.f2);
iDY.vec_ph2.f1 = min(iDY.intervalo.f1):(max(iDY.intervalo.f1)-min(iDY.intervalo.f1))/pr2:max(iDY.intervalo.f1);

iDY.fPo.f1 = matlabFunction(iDY.iME.f1*Vo, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
iDY.fPo.f2 = matlabFunction(iDY.iME.f2*Vo, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

figure
hold on
plot(iDY.vec_ph2.f1*180/pi,iDY.fPo.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,1,fs_num,iDY.vec_ph2.f1),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.vec_ph2.f1*180/pi,iDY.fPo.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,1/2,fs_num,iDY.vec_ph2.f1),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(iDY.vec_ph2.f2*180/pi,iDY.fPo.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,1,fs_num,iDY.vec_ph2.f2),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(iDY.vec_ph2.f2*180/pi,iDY.fPo.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,1/2,fs_num,iDY.vec_ph2.f2),'Color',jetcustom(2,:),'LineWidth',1.5)
yline(pt_po)

legend({'400','300'},'Location','southeast','FontSize', 16)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
title('iDY')


%% potencia de saida em funcao do angulo
pr2 = pr*100;
YY.vec_ph2.f2 = min(YY.intervalo.f2):(max(YY.intervalo.f2)-min(YY.intervalo.f2))/pr2:max(YY.intervalo.f2);
YY.vec_ph2.f1 = min(YY.intervalo.f1):(max(YY.intervalo.f1)-min(YY.intervalo.f1))/pr2:max(YY.intervalo.f1);

YY.fPo.f1 = matlabFunction(YY.iME.f1*Vo, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fPo.f2 = matlabFunction(YY.iME.f2*Vo, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

figure
hold on
plot(YY.vec_ph2.f1*180/pi,YY.fPo.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,1,fs_num,YY.vec_ph2.f1),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.vec_ph2.f1*180/pi,YY.fPo.f1(L1_num,L2_num,Ldab_num,M_num,Vi_num,3/4,fs_num,YY.vec_ph2.f1),'Color',jetcustom(2,:),'LineWidth',1.5)

plot(YY.vec_ph2.f2*180/pi,YY.fPo.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,1,fs_num,YY.vec_ph2.f2),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(YY.vec_ph2.f2*180/pi,YY.fPo.f2(L1_num,L2_num,Ldab_num,M_num,Vi_num,3/4,fs_num,YY.vec_ph2.f2),'Color',jetcustom(2,:),'LineWidth',1.5)
yline(pt_po)
legend({'400','300'},'Location','southeast','FontSize', 16)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
title('YY')


% 
% 
% %% primeiro trafo DiY
% load('YY.mat');
% load('YD.mat');
% load('DfD.mat');
% load('DiD.mat');
% load('DiY.mat');
% load('DfY.mat');
% 
% % Combine the structs into a single struct
% l.eq = struct('YY', YY, 'YD', YD, 'DfD', DfD, ...
%               'DiD', DiD, 'DiY', DiY, 'DfY', DfY);
% 
% trafoo = "DfY";
% %% os plots, ponto A 
% l.pr.phi = deg2rad(-15);
% [out] = f_equations_plot(l,trafoo);
% C = num2cell(out);
% vPP(1:3,:) = cell2mat(C(1:3,:));
% vSS(1:3,:) = cell2mat(C(4:6,:));
% vLLdab(1:3,:) = cell2mat(C(7:9,:));
% hbb(1:3,:) = cell2mat(C(10:12,:));
% HBB(1:3,:) = cell2mat(C(13:15,:));
% idab_plot(1:3,:) = cell2mat(C(16:18,:));
% x0s(1:6,:) = cell2mat(C(19:24,:));
% ts(1,:) = cell2mat(C(25,:));
% sec_switch_plot(1,:) = cell2mat(C(26,:));
% sf_plot(1:6,:) = cell2mat(C(27:32,:));




