clear; close all; clc;

syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po pi dt t real

phi_num = deg2rad(-90);  %[MUDAR]


[f_IL_rms,f_Itrf_sec_rms,f_Iin_med,f_Iin_rms,f_Iout_med,f_Iout_rms,f_phi,f_Isw_p_rms,f_Isw_s_rms,f_Ip,f_Is,f_IL_rms_c_k,f_Itrf_sec_rms_c_k] = dabDinY(phi_num);


Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
Ld2_num = 2e-6;
Lm_num = 100e-6;
n_num = 5/9;
Vo_num = d_num*Vi_num;
Po_num = 3000;


f_IL_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Itrf_sec_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_IL_rms_c_k(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num, 7)*sqrt(2)
f_Itrf_sec_rms_c_k(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num, 7)*sqrt(2)
f_Iin_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iin_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iout_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iout_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_phi(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), Po_num, fs_num, Vi_num, Vo_num)*57.295779513082323
f_Isw_p_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Isw_s_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Ip(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Is(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)

