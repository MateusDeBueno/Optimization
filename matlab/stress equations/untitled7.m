clear; close all; clc;

syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po pi dt t real

phi_num = deg2rad(-80);  %[MUDAR]

file = '0to60';

[f_IL_rms,f_Itrf_sec_rms,f_Iin_med,f_Iin_rms,f_Iout_med,f_Iout_rms,f_phi,f_Isw_p_rms,f_Isw_s_rms,f_Ip,f_Is,f_IL_rms_c_k,f_Itrf_sec_rms_c_k] = dabYY2(phi_num);


Vi = 400;
d = 1;
fs = 100e3;
Ldab = 61e-6;
Ld1 = 2e-6;
Ld2 = 2e-6;
Lm = 700e-6;
phi = deg2rad(-80);  %[MUDAR]
n = 5/9;
Vo = d*Vi;
Po = 3000;
pi = 3.141592653589793;