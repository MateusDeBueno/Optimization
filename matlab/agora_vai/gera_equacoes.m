clear; close all; clc;

addpath('utils')
addpath('utils_transf')
addpath('utils_loss')

syms Ld1 Ld2 n Lm Po t L1 L2 Ldab M fs Vi d dt a b ki N Ac Ve real positive
syms phi real

% parpool('myProf',8)

phi_num = deg2rad(-150);  %[MUDAR]
[YY.eq1.f,YY.eq1.fh] = full_YY(phi_num);
phi_num = deg2rad(-90);  %[MUDAR]
[YY.eq2.f,YY.eq2.fh] = full_YY(phi_num);
phi_num = deg2rad(-30);  %[MUDAR]
[YY.eq3.f,YY.eq3.fh] = full_YY(phi_num);
phi_num = deg2rad(30);  %[MUDAR]
[YY.eq4.f,YY.eq4.fh] = full_YY(phi_num);
phi_num = deg2rad(90);  %[MUDAR]
[YY.eq5.f,YY.eq5.fh] = full_YY(phi_num);
phi_num = deg2rad(150);  %[MUDAR]
[YY.eq6.f,YY.eq6.fh] = full_YY(phi_num);
save('YY.mat', 'YY')

phi_num = deg2rad(-150);  %[MUDAR]
[DiY.eq1.f,DiY.eq1.fh] = full_DiY(phi_num);
phi_num = deg2rad(-90);  %[MUDAR]
[DiY.eq2.f,DiY.eq2.fh] = full_DiY(phi_num);
phi_num = deg2rad(-30);  %[MUDAR]
[DiY.eq3.f,DiY.eq3.fh] = full_DiY(phi_num);
phi_num = deg2rad(30);  %[MUDAR]
[DiY.eq4.f,DiY.eq4.fh] = full_DiY(phi_num);
phi_num = deg2rad(90);  %[MUDAR]
[DiY.eq5.f,DiY.eq5.fh] = full_DiY(phi_num);
phi_num = deg2rad(150);  %[MUDAR]
[DiY.eq6.f,DiY.eq6.fh] = full_DiY(phi_num);
save('DiY.mat', 'DiY')

phi_num = deg2rad(-150);  %[MUDAR]
[DfY.eq1.f,DfY.eq1.fh] = full_DfY(phi_num);
phi_num = deg2rad(-90);  %[MUDAR]
[DfY.eq2.f,DfY.eq2.fh] = full_DfY(phi_num);
phi_num = deg2rad(-30);  %[MUDAR]
[DfY.eq3.f,DfY.eq3.fh] = full_DfY(phi_num);
phi_num = deg2rad(30);  %[MUDAR]
[DfY.eq4.f,DfY.eq4.fh] = full_DfY(phi_num);
phi_num = deg2rad(90);  %[MUDAR]
[DfY.eq5.f,DfY.eq5.fh] = full_DfY(phi_num);
phi_num = deg2rad(150);  %[MUDAR]
[DfY.eq6.f,DfY.eq6.fh] = full_DfY(phi_num);
save('DfY.mat', 'DfY')

phi_num = deg2rad(-150);  %[MUDAR]
[DfD.eq1.f,DfD.eq1.fh] = full_DfD(phi_num);
phi_num = deg2rad(-90);  %[MUDAR]
[DfD.eq2.f,DfD.eq2.fh] = full_DfD(phi_num);
phi_num = deg2rad(-30);  %[MUDAR]
[DfD.eq3.f,DfD.eq3.fh] = full_DfD(phi_num);
phi_num = deg2rad(30);  %[MUDAR]
[DfD.eq4.f,DfD.eq4.fh] = full_DfD(phi_num);
phi_num = deg2rad(90);  %[MUDAR]
[DfD.eq5.f,DfD.eq5.fh] = full_DfD(phi_num);
phi_num = deg2rad(150);  %[MUDAR]
[DfD.eq6.f,DfD.eq6.fh] = full_DfD(phi_num);
save('DfD.mat', 'DfD')

phi_num = deg2rad(-150);  %[MUDAR]
[DiD.eq1.f,DiD.eq1.fh] = full_DiD(phi_num);
phi_num = deg2rad(-90);  %[MUDAR]
[DiD.eq2.f,DiD.eq2.fh] = full_DiD(phi_num);
phi_num = deg2rad(-30);  %[MUDAR]
[DiD.eq3.f,DiD.eq3.fh] = full_DiD(phi_num);
phi_num = deg2rad(30);  %[MUDAR]
[DiD.eq4.f,DiD.eq4.fh] = full_DiD(phi_num);
phi_num = deg2rad(90);  %[MUDAR]
[DiD.eq5.f,DiD.eq5.fh] = full_DiD(phi_num);
phi_num = deg2rad(150);  %[MUDAR]
[DiD.eq6.f,DiD.eq6.fh] = full_DiD(phi_num);
save('DiD.mat', 'DiD')

phi_num = deg2rad(-150);  %[MUDAR]
[YD.eq1.f,YD.eq1.fh] = full_YD(phi_num);
phi_num = deg2rad(-90);  %[MUDAR]
[YD.eq2.f,YD.eq2.fh] = full_YD(phi_num);
phi_num = deg2rad(-30);  %[MUDAR]
[YD.eq3.f,YD.eq3.fh] = full_YD(phi_num);
phi_num = deg2rad(30);  %[MUDAR]
[YD.eq4.f,YD.eq4.fh] = full_YD(phi_num);
phi_num = deg2rad(90);  %[MUDAR]
[YD.eq5.f,YD.eq5.fh] = full_YD(phi_num);
phi_num = deg2rad(150);  %[MUDAR]
[YD.eq6.f,YD.eq6.fh] = full_YD(phi_num);
save('YD.mat', 'YD')