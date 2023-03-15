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

l.off = load('f_fitted_off.mat');
l.on = load('f_fitted_on.mat');

color1 = [0.045 0.245 0.745]; % blue
color2 = [0.635 0.635 0.635]; % gray

syms Ld1 Ld2 n Lm Po t L1 L2 Ldab M fs Vi d dt a b ki N Ac Ve real positive
syms phi real

Vo = Vi*d;
Ts = 1/fs;

phi_num = deg2rad(50);  %[MUDAR]
l.eq1.output = full_YY(phi_num);
phi_num = deg2rad(70);  %[MUDAR]
l.eq2.output = full_YY(phi_num);


l.pr.phi = deg2rad(67);
l.pr.Vi = 400;
l.pr.d = 1;
l.pr.fs = 100e3;
l.pr.Ldab = 61e-6;
l.pr.Ld1 = 2e-6;
l.pr.n = 5/9;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 700e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);

l.L.a = 1.394;
l.L.b = 2.248;
l.L.kc = 2.448;
l.L.int_ki = integral(@(theta) abs(cos(theta)).^l.L.a,0,2*pi);
l.L.ki = l.L.kc/(2^(l.L.b-l.L.a)*(2*pi)^(l.L.a-1)*l.L.int_ki);
l.L.Ac = 421.3e-6;
l.L.Ve = 52.1e-6;
l.L.N = 14;

l.tr.a = 1.585;
l.tr.b = 2.756;
l.tr.kc = 0.402;
l.tr.int_ki = integral(@(theta) abs(cos(theta)).^l.tr.a,0,2*pi);
l.tr.ki = l.tr.kc/(2^(l.tr.b-l.tr.a)*(2*pi)^(l.tr.a-1)*l.tr.int_ki);
l.tr.Ac = 683e-6;
l.tr.Ve = 146.79e-6;
l.tr.N = 9;

x = l.eq2.output(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi, ...
    l.L.ki,l.L.b,l.L.a,l.L.Ac,l.L.N,l.L.Ve, ...
    l.tr.ki,l.tr.b,l.tr.a,l.tr.Ac,l.tr.N,l.tr.Ve);

C = num2cell(x);
[hbrm,HBrm,Ip,Is,Pm,idrm,ilrm,iLrm,iSwPrm,iSwSrm,P_core_tr,Bpk_tr,P_core_L,Bpk_L] = C{:};


%%

fts = matlabFunction(ts, 'vars', {fs,phi});
fiSwPrm = matlabFunction(iSwPrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwSrm = matlabFunction(iSwSrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIp = matlabFunction(Ip, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs = matlabFunction(Is, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fP_tr = matlabFunction(P_tr, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi,a,b,ki,N,Ac,Ae});
fP_L = matlabFunction(P_L, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi,a,b,ki,N,Ac,Ae});

%%
iSwPrm_num = fiSwPrm(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi_num);
iSwSrm_num = fiSwSrm(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi_num);
Ip_num = fIp(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi_num);
Is_num = fIs(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi_num);


%loss calculation
[cond_loss] = f_cond_loss(iSwPrm_num,iSwSrm_num);
prim_loss = f_sw_loss(fs_num, Vi_num, Ip_num, l)*6;
sec_loss = f_sw_loss(fs_num, Vi_num*d_num, Is_num, l)*6;
p_core_L = 3*fPv_L(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi_num,l.L.a,l.L.b,l.L.ki,l.L.N,l.L.Ac)*l.L.Ve;
p_core_tr = fPv_tr(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi_num,l.tr.a,l.tr.b,l.tr.ki,l.tr.N,l.tr.Ac)*l.tr.Ve;
