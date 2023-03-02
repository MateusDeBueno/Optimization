
clear; close all; clc;

addpath('utils')
addpath('utils_transf')

syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po pi k dt t real

%% SALVAR EQUACOES DO DinY
ii = 0;
for phi_num=-150:60:150
    phi_num
    ii = ii + 1;
    [str(ii,:),str_h(ii,:),fphi(ii,:)] = dabDinY(deg2rad(phi_num));
end

trafo.DinY.str.f1 = matlabFunction(str(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinY.str.f2 = matlabFunction(str(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinY.str.f3 = matlabFunction(str(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinY.str.f4 = matlabFunction(str(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinY.str.f5 = matlabFunction(str(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinY.str.f6 = matlabFunction(str(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

trafo.DinY.str_h.f1 = matlabFunction(str_h(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinY.str_h.f2 = matlabFunction(str_h(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinY.str_h.f3 = matlabFunction(str_h(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinY.str_h.f4 = matlabFunction(str_h(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinY.str_h.f5 = matlabFunction(str_h(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinY.str_h.f6 = matlabFunction(str_h(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});

trafo.DinY.fphi.f1 = matlabFunction(fphi(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});
trafo.DinY.fphi.f2 = matlabFunction(fphi(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});
trafo.DinY.fphi.f3 = matlabFunction(fphi(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});
trafo.DinY.fphi.f4 = matlabFunction(fphi(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});
trafo.DinY.fphi.f5 = matlabFunction(fphi(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});
trafo.DinY.fphi.f6 = matlabFunction(fphi(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});





Po_num = -3700;


opc1 = trafo.DinY.fphi.f1(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)*180/3.1415
opc2 = trafo.DinY.fphi.f2(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)*180/3.1415
opc3 = trafo.DinY.fphi.f3(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)*180/3.1415
opc4 = trafo.DinY.fphi.f4(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)*180/3.1415
opc5 = trafo.DinY.fphi.f5(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)*180/3.1415
opc6 = trafo.DinY.fphi.f6(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)*180/3.1415








%% SALVAR EQUACOES DO DinY
ii = 0;
for phi_num=-150:60:150
    phi_num
    ii = ii + 1;
    [str(ii,:),str_h(ii,:),fphi(ii,:)] = dabYY(deg2rad(phi_num));
end

trafo.YY.str.f1 = matlabFunction(str(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YY.str.f2 = matlabFunction(str(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YY.str.f3 = matlabFunction(str(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YY.str.f4 = matlabFunction(str(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YY.str.f5 = matlabFunction(str(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YY.str.f6 = matlabFunction(str(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

trafo.YY.str_h.f1 = matlabFunction(str_h(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YY.str_h.f2 = matlabFunction(str_h(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YY.str_h.f3 = matlabFunction(str_h(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YY.str_h.f4 = matlabFunction(str_h(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YY.str_h.f5 = matlabFunction(str_h(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YY.str_h.f6 = matlabFunction(str_h(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});

%% SALVAR EQUACOES DO DinY
ii = 0;
for phi_num=-150:60:150
    phi_num
    ii = ii + 1;
    [str(ii,:),str_h(ii,:),fphi(ii,:)] = dabYD(deg2rad(phi_num));
end

trafo.YD.str.f1 = matlabFunction(str(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YD.str.f2 = matlabFunction(str(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YD.str.f3 = matlabFunction(str(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YD.str.f4 = matlabFunction(str(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YD.str.f5 = matlabFunction(str(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.YD.str.f6 = matlabFunction(str(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

trafo.YD.str_h.f1 = matlabFunction(str_h(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YD.str_h.f2 = matlabFunction(str_h(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YD.str_h.f3 = matlabFunction(str_h(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YD.str_h.f4 = matlabFunction(str_h(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YD.str_h.f5 = matlabFunction(str_h(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.YD.str_h.f6 = matlabFunction(str_h(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});

%% SALVAR EQUACOES DO DinY
ii = 0;
for phi_num=-150:60:150
    phi_num
    ii = ii + 1;
    [str(ii,:),str_h(ii,:),fphi(ii,:)] = dabDinD(deg2rad(phi_num));
end

trafo.DinD.str.f1 = matlabFunction(str(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinD.str.f2 = matlabFunction(str(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinD.str.f3 = matlabFunction(str(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinD.str.f4 = matlabFunction(str(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinD.str.f5 = matlabFunction(str(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
trafo.DinD.str.f6 = matlabFunction(str(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

trafo.DinD.str_h.f1 = matlabFunction(str_h(1,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinD.str_h.f2 = matlabFunction(str_h(2,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinD.str_h.f3 = matlabFunction(str_h(3,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinD.str_h.f4 = matlabFunction(str_h(4,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinD.str_h.f5 = matlabFunction(str_h(5,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
trafo.DinD.str_h.f6 = matlabFunction(str_h(6,:), 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});

%%

save('trafo.mat','trafo')

%%
clear; close all; clc;
load('trafo.mat')

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 10e-6;
Ld2_num = 5e-6;
Lm_num = 700e-6;
phi_num = deg2rad(-150);  %[MUDAR]
n_num = 5/9;
Vo_num = d_num*Vi_num;
k_num = 1;
Po_num = -4700;


trafo.DinY.fphi.f1(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)*180/3.14

[out] = fDinY(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num)
[out] = fDinY_harm(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num)
[out] = fYY(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num)
[out] = fYY_harm(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num)