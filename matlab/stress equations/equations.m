clear; close all; clc;

% syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po pi dt t real


phi_num = deg2rad(50)
output = dabDinY(phi_num)



syms x
f = @(x) (x>=0)*x.^2 + (x<0)*x
f1 = @(x) (x>=0)*x.^10 + (x<0)*5*x

vec = {f, f1}
valueOfSecondFunction = vec{1}(0:2)


%%




load("YYoutput.mat")

f1 = output{1}(1)

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


output{1}(1)
(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)





%%
% count = 1;
% for phi_num=-150:60:150
%     output{count,:} = dabDinY(deg2rad(phi_num));
%     count = count + 1;
% end
% save("YYoutput.mat","output")


% 
%     if (-Pi <= phi_num && phi_num <= -2*Pi/3)
%         situacao = 0;
%     elseif (-2*Pi/3 < phi_num && phi_num <= -Pi/3)
% 
% f = @(Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo) (-pi<=phi && phi<=-2*pi/3)*output{1}(1) + (-2*pi/3<=phi && phi<=-pi/3)*output{2}(1)
% 
% 
% %%
% 
% f = @(phi) (0<=phi && phi<=pi/6)*output{1} + (phi<0)*output{2}
% 
% 
% %%
% 
% f = @(Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo) (0<=phi && phi<=pi/6)*f0to60.f_IL_rms + (phi<0)*f0to60.f_IL_rms
% 
% 
% %%
% 
% save('f0to60.mat','f0to60')
% 
% load('f0to60.mat')
% 
% phi_num = deg2rad(-90);  %[MUDAR]
% % [f_IL_rms,f_Itrf_sec_rms,f_Iin_med,f_Iin_rms,f_Iout_med,f_Iout_rms,f_phi,f_Isw_p_rms,f_Isw_s_rms,f_Ip,f_Is,f_IL_rms_c_k,f_Itrf_sec_rms_c_k] = dabDinY(phi_num);
% 
% % f = @(x) (x>=0)*x.^2 + (x<0)*x
% 
% Vi_num = 400;
% d_num = 1;
% fs_num = 100e3;
% Ldab_num = 61e-6;
% Ld1_num = 2e-6;
% Ld2_num = 2e-6;
% Lm_num = 100e-6;
% n_num = 5/9;
% Vo_num = d_num*Vi_num;
% Po_num = 3000;
% 
% f(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% 
% 
% f_IL_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Itrf_sec_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_IL_rms_c_k(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num, 7)*sqrt(2)
% f_Itrf_sec_rms_c_k(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num, 7)*sqrt(2)
% f_Iin_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Iin_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Iout_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Iout_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_phi(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), Po_num, fs_num, Vi_num, Vo_num)*57.295779513082323
% f_Isw_p_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Isw_s_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Ip(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Is(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)

