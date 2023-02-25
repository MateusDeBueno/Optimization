clear; close all; clc;

syms phi fs Vi Vo real
syms Ldab Ld1 Ld2 Lm n Po real
syms Po real

addpath('utils')

% Parametros para validacao numerica
Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
Ld2_num = 2e-6;
Lm_num = [700e-6, 1e-3, 1000e-3];
phi_num = deg2rad(80);
n_num = 5/9;
Vo_num = d_num*Vi_num;
Po_num = 3500;
variaveis =         [Ldab,      Ld1,        Ld2,        Lm,         Vi,         Vo,         fs,         n,          phi];
ponto_de_operacao = [Ldab_num,  Ld1_num,    Ld2_num,    Lm_num(1),  Vi_num,     Vo_num,     fs_num,     n_num,      phi_num];


%Transformada de clark
Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1];
Tclm = Tcl(1:2,:); %ignorar nivel zero
Tclx = kron(eye(2), Tclm);

Ts = 1/fs;

%Definicao do angulo
angs = [0 phi-pi/3 pi/3 phi 2*pi/3 pi/3+phi];
angs = [angs,pi+angs,2*pi];
ts = angs*Ts/(2*pi);

%Definicao dos estados de comutacao, ponto medio no meio do barramento
sf = -0.5+[1 0 1; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]';
sf = kron(sf, [1, 1]); %abc
sf = [sf; sf(:,end-2:end) sf(:,1:end-3)]; %a1b1c1 to a2b2c2
scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos

%%  Validação dos Pulsos de Comutação
figure(1)
subplot(2,1,1)
stairs(tf(phi_num, fs_num) * 1e6, ([sf sf(:,1)]'+[0 -2 -4 -6 -8 -10]), 'LineWidth', 2)

xlim([0 1/fs_num]*1e6)
ylabel('s')
xlabel('t [us]')
legend('Primary', 'Secondary')
grid on
subplot(2,1,2)
stairs(tf(phi_num, fs_num) * 1e6, ([scl scl(:,1)]'+[0 -2 -4 -6]), 'LineWidth', 2)
xlim([0 1/fs_num]*1e6)
ylabel('s_ab')
xlabel('t [us]')
legend('Primary', 'Secondary')
grid on

%% Resolvendo circuito monofasico equivalente
L = [Ldab+Ld1; Ld2/n^2];
Udc = [Vi; Vo];
u = Udc.*[1; 1/n];

x = sym('x_', [2,1], 'real');
dx = sym('dx_', [2,1], 'real'); %dx1 -> iLdab ou iLl1 (em serie) dx2 -> iLl2
eq = sym('eq', [2,1], 'real');
s = sym('s', [2,1], 'real');

eq(1) = u(1)*s(1) == L(1)*dx(1) + Lm*(dx(1)-dx(2)*n);
eq(2) = u(2)*s(2) == Lm*(dx(1)-dx(2)*n) - L(2)*dx(2)*n;

dxs = simplify(struct2array(solve(eq, dx))).';
M = simplify(equationsToMatrix(dxs, [x; s]));
A = M(:,1:2);
B = M(:,3:4); %devidas de iLdab e iLd

A = kron(A, eye(2));
B = kron(B, eye(2));

%% Obter valores de regime permanente

x = sym('x_', [4,1], 'real');
syms dt real

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

%% Corrente eficaz nos estados
syms t real

derivadas = pinv(Tclx)*B*scl;

[x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,(1:12),(1:12));

IL_rms = x_rms(1);
Itrf_sec_rms = x_rms(4);

f_IL_rms = matlabFunction(IL_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Itrf_sec_rms = matlabFunction(Itrf_sec_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
    
%% Corrente de entrada

[x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,(3:4),(3:4));
Iin_med = x_mean(1);
Iin_rms = x_rms(1);

f_Iin_med = matlabFunction(Iin_med, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Iin_rms = matlabFunction(Iin_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% Corrente de saida

[x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,(6:7),(6:7));
Iout_med = x_mean(4);
Iout_rms = x_rms(4);

f_Iout_med = matlabFunction(Iout_med, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Iout_rms = matlabFunction(Iout_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

power_equation = Iout_med==Po/Vo;
fphi = rhs(isolate(power_equation,phi));
phi_solutions = solve(power_equation,phi);

phi_solution = phi_solutions(2);
f_phi = matlabFunction(phi_solution, 'Vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi, Vo});

% f_phi(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), 3800, fs_num, Vi_num, Vo_num)*180/pi
%% Corrente nas chaves
[x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,(1:6),(1:12));

Isw_p_rms = x_rms(1);
Isw_s_rms = x_rms(4);

f_Isw_p_rms = matlabFunction(Isw_p_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Isw_s_rms = matlabFunction(Isw_s_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% Corrente de comutacao

%primario
Ip = x0s(1,1);

%secundario
Is = -x0s(4,4);

f_Is = matlabFunction(Is, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Ip = matlabFunction(Ip, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% valores

f_derivadas = matlabFunction(derivadas, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});


f_IL_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)


%% fourier currents tentando plotar
% % % % derivadass = f_derivadas(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num);
% % % % x0ss = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num);
% % % % time_intervals = tf(phi_num,fs_num);
% % % % time = 0:1/(fs_num*1000):1/fs_num;
% % % % count = 0;
% % % % for t = time
% % % %     count = count + 1;
% % % %     
% % % %     for i=1:length(time_intervals)-1
% % % %         if (time_intervals(i) <= t && t < time_intervals(i+1))
% % % %             interval = i;
% % % %         end
% % % %     end
% % % %     t_eq = t - time_intervals(interval);
% % % %     Ia(count) = derivadass(1,interval)*t_eq + x0ss(interval,1);
% % % % end
% % % % 
% % % % 
% % % % 
% % % % figure
% % % % plot(time,Ia)
% % % % 
% % % % for i=1:length(ts)-1
% % % %     dt = ts(i+1) - ts(i);
% % % %     equation_intervalo(:,i) = (derivadas(:,i).*ones(6,1)*t+x0s(:,i));
% % % % end

%%


% syms t k real
% 
% a_k = 0;
% for i=1:length(ts)-1
%     dt = ts(i+1) - ts(i);
%     equation_intervalo = (derivadas(:,i).*ones(6,1)*t+x0s(:,i));
%     a_k = a_k + int(equation_intervalo*cos(2*pi*k*t/Ts), [0, dt]);
% end
% a_k = simplify(a_k*2/Ts);
% 
% b_k = 0;
% for i=1:length(ts)-1
%     dt = ts(i+1) - ts(i);
%     equation_intervalo = (derivadas(:,i).*ones(6,1)*t+x0s(:,i));
%     b_k = b_k + int(equation_intervalo*sin(2*pi*k*t/Ts), [0, dt]);
% end
% b_k = simplify(b_k*2/Ts);
% 
% c_k = abs(a_k-1i*b_k); %cos e sin to exponential
% c_k = simplify(c_k);
% c_k_rms = c_k/sqrt(2); %peak to rms
% 
% 
% f_c_k = matlabFunction(c_k', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});





%% equacoes de saidas
output = [IL_rms,Itrf_sec_rms,Iin_med,Iin_rms,Iout_med,Iout_rms,fphi,Isw_p_rms,Isw_s_rms,Ip,Is]';
vpa(subs(output,variaveis,ponto_de_operacao))

% f_IL_rms
% f_Itrf_sec_rms
% f_Iin_med
% f_Iin_rms
% f_Iout_med
% f_Iout_rms
% f_phi
% f_Isw_p_rms
% f_Isw_s_rms
% f_Is
% f_Ip

% f_YY_maior60 = matlabFunction(M, 'vars', variaveis);
save('dabYY_functions_maior60.mat', "f_IL_rms", "f_Itrf_sec_rms", "f_Iin_med", "f_Iin_rms", ...
    "f_Iout_med", "f_Iout_rms", "f_phi", "f_Isw_p_rms", "f_Isw_s_rms", "f_Is", "f_Ip", '-mat')


%% Plots

figure
hold on
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'r','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'b','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'g')
hold off
grid on
grid minor
xlim([0 1/fs_num])


%% Plot corrente secundario
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

%% Plot corrente da fase A, entrada e saida
M_fase_a = [1 0; 0 0; 0 0; 0 1; 0 0; 0 0];

figure
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*M_fase_a,'LineWidth',2)
grid on
grid minor
xlim([0 1/fs_num])
