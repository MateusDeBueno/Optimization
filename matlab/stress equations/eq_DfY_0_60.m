clear; close all; clc;
syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po dt t real

% Parametros para validacao numerica
Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
Ld2_num = 2e-6;
Lm_num = [100e-6, 1e-3, 1000e-3];
phi_num = deg2rad(-90);  %[MUDAR]
n_num = 5/9;
Vo_num = d_num*Vi_num;
Po_num = 3000;
pi_num = 3.141592653589793;



syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n Po dt k t real
    
pi_num = 3.141592653589793;
    
Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
Tclm = Tcl(1:2,:); %ignorar nivel zero
Tclx = kron(eye(2), Tclm);
    
Ts = 1/fs;

[sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi);

sf_line(1,:) = sf(1,:)-sf(2,:);
sf_line(2,:) = sf(2,:)-sf(3,:);
sf_line(3,:) = sf(3,:)-sf(1,:);

sf = [sf_line; sf_s]; %a1b1c1 to a2b2c2

scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
ts = simplify(ang)*Ts/(2*pi);
tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos
        
%% Resolvendo circuito monofasico equivalente [MUDAR]
L = [Ldab+Ld1/3; Ld2/n^2];
Udc = [Vi; Vo];
u = Udc.*[1; 1/n];

x = sym('x_', [2,1], 'real');
dx = sym('dx_', [2,1], 'real'); %dx1 -> iLdab ou iLl1 (em serie) dx2 -> iLl2
eq = sym('eq', [2,1], 'real');
s = sym('s', [2,1], 'real');

eq(1) = u(1)*s(1) == L(1)*dx(1) + Lm/sqrt(3)*(dx(1)-dx(2)*n);
eq(2) = u(2)*s(2) == Lm/sqrt(3)*(dx(1)-dx(2)*n) - L(2)*dx(2)*n;

dxs = simplify(struct2array(solve(eq, dx))).';
M = equationsToMatrix(dxs, [x; s]);
A = M(:,1:2);
B = M(:,3:4); %derivadas de iLdab e iLd

A = kron(A, eye(2));
B = kron(B, eye(2));
    
%% Obter valores de regime permanente
g = simplify(expm(A*dt));
h = B*dt;
x = sym('x_', [4,length(ts)], 'real');
x0 = sym('x0_',[4,1], 'real');

x(:,1) = x0;
for i=1:length(ts)-1
    dts = (ts(i+1)-ts(i));
    x(:,i+1) = simplify(subs(g, dt, dts)*x(:,i) + subs(h, dt, dts)*scl(:,i));
end

x0x = struct2array(solve(x(:,1) == -x(:,7), x0)).';
x0s = simplify(pinv(Tclx)*subs(x, x0, x0x));
iSw = matlabFunction(x0s.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
derivadas = pinv(Tclx)*B*scl; %derivadas dos estados, indutor e trafo secundario

%% Corrente eficaz nos estados
[x_rms,x_mean] = rms_and_mean(derivadas,x0s,ts,(1:12),(1:12));
IL_rms = x_rms(1);
Itrf_sec_rms = x_rms(4);
f_IL_rms = matlabFunction(IL_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Itrf_sec_rms = matlabFunction(Itrf_sec_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% Fourier dos estados
% IL_rms_c_k = ck_fourier(ts,derivadas(1,:),x0s(1,:));
% Itrf_sec_rms_c_k = ck_fourier(ts,derivadas(4,:),x0s(4,:));
% 
% f_IL_rms_c_k = matlabFunction(IL_rms_c_k, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});
% f_Itrf_sec_rms_c_k = matlabFunction(Itrf_sec_rms_c_k, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo, k});

%% Corrente do half bridge do primario [MUDAR]
derivadas_hb_p = derivadas([1 2 3],:) - derivadas([3 1 2],:);
pts_inics_hb_p = x0s([1 2 3],:) - x0s([3 1 2],:);

%% Corrente do half bridge do primario [MUDAR]
derivadas_hb_s = derivadas(4:6,:);
pts_inics_hb_s = x0s(4:6,:);

%% Corrente de entrada
target = [1;0;0];
[etapas] = pega_etapa(sf_p,target);
[x_rms,x_mean] = rms_and_mean(derivadas_hb_p(1,:),pts_inics_hb_p(1,:),ts,etapas,etapas);
f_Iin_med = matlabFunction(x_mean, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Iin_rms = matlabFunction(x_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% Corrente de saida
target = [1;0;0];
[etapas] = pega_etapa(sf_s,target);
[x_rms,x_mean] = rms_and_mean(derivadas_hb_s(1,:),pts_inics_hb_s(1,:),ts,etapas,etapas);
f_Iout_med = matlabFunction(x_mean, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Iout_rms = matlabFunction(x_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

power_equation = x_mean==Po/Vo;
phi_solutions = solve(power_equation,phi);
f_phi = matlabFunction(phi_solutions, 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});

%% Corrente nas chaves
[x_rms,x_mean] = rms_and_mean(derivadas_hb_p(1,:),pts_inics_hb_p(1,:),ts,(1:6),(1:12));
Isw_p_rms = x_rms;

[x_rms,x_mean] = rms_and_mean(derivadas_hb_s(1,:),pts_inics_hb_s(1,:),ts,(1:6),(1:12));
Isw_s_rms = x_rms;

f_Isw_p_rms = matlabFunction(Isw_p_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Isw_s_rms = matlabFunction(Isw_s_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% corrente de comutacao
%primario
Ip = pts_inics_hb_p(1,1);
f_Ip = matlabFunction(Ip, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%secundario
Is = -pts_inics_hb_s(1,sec_switch);
f_Is = matlabFunction(Is, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% plots
figure
hold on
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'r','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'b','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'g')
hold off
grid on
grid minor
xlim([0 1/fs_num])


%% plot corrente secundario
M_sec = [zeros(3,3),eye(3,3)]';
figure
hold on
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*M_sec,'r','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*M_sec,'b','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*M_sec,'g')
hold off
grid on
grid minor
xlim([0 1/fs_num])

%% plot corrente da fase A, entrada e saida
M_fase_a = [1 0; 0 0; 0 0; 0 1; 0 0; 0 0];
figure
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*M_fase_a,'LineWidth',2)
grid on
grid minor
xlim([0 1/fs_num])





f_IL_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Itrf_sec_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_IL_rms_c_k(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num, 1)
% f_Itrf_sec_rms_c_k(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num, 1)
f_Iin_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iin_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iout_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iout_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_phi(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), Po_num, fs_num, Vi_num, Vo_num)*57.295779513082323
f_Isw_p_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Isw_s_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Ip(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Is(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
