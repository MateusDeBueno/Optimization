clear; close all; clc;

syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po dt k t real


phi_num = deg2rad(-80);  %[MUDAR]
[IL_rms,Itrf_sec_rms,Iin_med,Iin_rms,Iout_med,Iout_rms,phi,Isw_p_rms,Isw_s_rms,Ip,Is,IL_rms_c_k,Itrf_sec_rms_c_k] = dabYY2(phi_num);

% func = matlabFunction(IL_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
save('equationss.mat', 'IL_rms','Itrf_sec_rms','Iin_med','Iin_rms','Iout_med','Iout_rms','phi','Isw_p_rms','Isw_s_rms','Ip','Is','IL_rms_c_k','Itrf_sec_rms_c_k')