clear; close all; clc;

addpath('utils')
addpath('utils_transf')

syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po pi k dt t real

phi_num = deg2rad(50);
[output,output_harmonic,output_phi] = dabDinD(phi_num)


Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
Ld2_num = 2e-6;
Lm_num = 700e-6;
n_num = 5/9;
Vo_num = d_num*Vi_num;
k_num = 1;


testa = matlabFunction(output, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});

% testa_phi(:,1) = matlabFunction(output_phi, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});



testa(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num)
