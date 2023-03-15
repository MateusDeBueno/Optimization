clear; close all; clc;

addpath('utils')
addpath('utils_transf')
addpath('utils_loss')

syms Ld1 Ld2 n Lm Po t L1 L2 Ldab M fs Vi d dt a b ki N Ac Ve real positive
syms phi real


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
