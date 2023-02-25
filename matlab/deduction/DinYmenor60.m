clear; close all; clc;

syms phi fs Vi Vo real
syms Ldab Ld1 Ld2 Lm n Po real
syms Po real

% Parametros para validacao numerica
Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
Ld2_num = 2e-6;
Lm_num = [700e-6, 1e-3, 1000e-3];
phi_num = deg2rad(50);
n_num = 5/9;
Vo_num = d_num*Vi_num;
Po_num = 3000;

variaveis =         [Ldab,      Ld1,        Ld2,        Lm,         Vi,         Vo,         fs,         n,          phi,        Po];
ponto_de_operacao = [Ldab_num,  Ld1_num,    Ld2_num,    Lm_num(1),  Vi_num,     Vo_num,     fs_num,     n_num,      phi_num,    Po_num];

%Transformada de clark
Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1];
Tclm = Tcl(1:2,:); %ignorar nivel zero
Tclx = kron(eye(2), Tclm);

Ts = 1/fs;

%Definicao do angulo
angs = [0 phi pi/3 pi/3+phi 2*pi/3 2*pi/3+phi];
angs = [angs,pi+angs,2*pi];
ts = angs*Ts/(2*pi);

%Definicao dos estados de comutacao, ponto medio no meio do barramento
sf = -0.5+[1 0 1; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]';
sf = kron(sf, [1, 1]); %abc


sf_line(1,:) = sf(1,:)-sf(2,:);
sf_line(2,:) = sf(2,:)-sf(3,:);
sf_line(3,:) = sf(3,:)-sf(1,:);


sf = [sf_line; sf(:, end) sf(:,1:end-1)]; %a1b1c1 to a2b2c2
scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos

%%  Validação dos Pulsos de Comutação
figure(1)
subplot(2,1,1)
stairs(tf(phi_num, fs_num) * 1e6, ([sf sf(:,1)]'+[0 -2 -4 -6 -8 -10]), 'LineWidth', 2)

xlim([0 1/fs_num]*1e6)
ylabel('s')
xlabel('t [\mus]')
legend('Primary', 'Secondary')
grid on
subplot(2,1,2)
stairs(tf(phi_num, fs_num) * 1e6, ([scl scl(:,1)]'+[0 -2 -4 -6]), 'LineWidth', 2)
xlim([0 1/fs_num]*1e6)
ylabel('s_{\alpha\beta}')
xlabel('t [\mus]')
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
M = equationsToMatrix(dxs, [x; s]);
A = M(:,1:2);
B = M(:,3:4); %devidas de iLdab e iLd

A = kron(A, eye(2));
B = kron(B, eye(2));

%% Obter valores de regime permanente

x0 = sym('x0_', [4,1], 'real');
x = sym('x_', [4,length(ts)], 'real');

x(:,1) = x0;
vs = [1 2 3 4 5 6 7 8 9 10 11 12];
for i=1:length(ts)-1
    dt = ts(i+1) - ts(i); %aqui tem uma simplificacao para somente retas
    x(:,i+1) = simplify(x(:,i) + dt*B*scl(:,vs(i))); %tempo x derivada x estados
end

x0x = struct2array(solve(x(:,1) == -x(:,7), x0)).'; %pega as 4 condicoes iniciais
x0x = subs(x, x0, x0x); %substitui as condicoes inicias nos estados

x0s = simplify(pinv(Tclx)*x0x);

iSw = matlabFunction(x0s.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% Corrente eficaz nos estados
syms t real

derivadas = pinv(Tclx)*B*scl;

[x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,(1:12),(1:12));

IL_rms = x_rms(1);
Itrf_sec_rms = x_rms(4);

f_IL_rms = matlabFunction(IL_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Itrf_sec_rms = matlabFunction(Itrf_sec_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});


%% Criar estado para corrente de saida do half bridge do primario

f_B = matlabFunction(B.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});


matriz_ponto_inicial = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num);
matriz_derivadas = f_B(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num);

[x_rms,x_mean] = f_rms_mean(simplify(derivadas(1,:)-derivadas(3,:)),simplify(x0s(1,:)-x0s(3,:)),ts,(3:4),(3:4));
Iin_med = x_mean(1)
Iin_rms = x_rms(1)

f_Iin_med = matlabFunction(Iin_med, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Iin_rms = matlabFunction(Iin_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

f_Iin_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iin_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)




figure
plot(tf(phi_num,fs_num), matriz_ponto_inicial(:,1)-matriz_ponto_inicial(:,3),'r','LineWidth',2)
grid on
grid minor
xlim([0 1/fs_num])



%% Corrente de entrada

[x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,(3:4),(3:4));
Iin_med = x_mean(1);
Iin_rms = x_rms(1);

f_Iin_med = matlabFunction(Iin_med, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Iin_rms = matlabFunction(Iin_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

f_Iin_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
f_Iin_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)



%% Corrente de saida

[x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,(4:5),(4:5));
Iout_med = x_mean(4);
Iout_rms = x_rms(4);


power_equation = Vo*Iout_med==Po;
f_Iout_med = matlabFunction(Iout_med, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Iout_rms = matlabFunction(Iout_rms, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

% f_Iout_med(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
% f_Iout_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)


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

%% corrente de comutacao

%primario
Ip = x0s(1,1);

%secundario
Is = -x0s(4,2);

f_Is = matlabFunction(Is, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
f_Ip = matlabFunction(Ip, 'Vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

%% valores

f_IL_rms(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)
%% fourier currents

% 
% syms t k real
% 
% a_k = 0;
% for i=1:length(ts)-1
%     dt = ts(i+1) - ts(i);
%     equation_intervalo = (derivadas(1,i)*t+x0s(1,i));
%     a_k = a_k + int(equation_intervalo*cos(2*pi*k*t/Ts), [ts(i), ts(i+1)]);
% end
% a_k = simplify(a_k*2/Ts);
% 
% b_k = 0;
% for i=1:length(ts)-1
%     equation_intervalo = (derivadas(1,i)*t+x0s(1,i));
%     b_k = b_k + int(equation_intervalo*sin(2*pi*k*t/Ts), [ts(i), ts(i+1)]);
% end
% b_k = simplify(b_k*2/Ts);
% 
% c_k = abs(a_k-1i*b_k); %cos e sin to exponential
% c_k = simplify(c_k);
% c_k_rms = c_k/sqrt(2); %peak to rms
% % 
% % 
% % f_c_k = matlabFunction(c_k', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
% 
% vpa(subs(subs(c_k,variaveis,ponto_de_operacao), k, 1))


%% equacoes de saidas
output = [IL_rms,Itrf_sec_rms,Iin_med,Iin_rms,Iout_med,Iout_rms,fphi,Isw_p_rms,Isw_s_rms,Ip,Is]';
vpa(subs(output,variaveis,ponto_de_operacao))

%% Lista de funcoes


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
save('dabYY_functions_menor60.mat', "f_IL_rms", "f_Itrf_sec_rms", "f_Iin_med", "f_Iin_rms", ...
    "f_Iout_med", "f_Iout_rms", "f_phi", "f_Isw_p_rms", "f_Isw_s_rms", "f_Is", "f_Ip", '-mat')


%% plots

figure
hold on
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'r','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(2), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'b','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(3), phi_num, fs_num, Vi_num, Vo_num)*eye(6,3),'g')
hold off
grid on
grid minor
xlim([0 1/fs_num])


%% plot corrente secundario
M_sec = [zeros(3,3),eye(3,3)]';


figure
hold on
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*M_sec,'r','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(2), phi_num, fs_num, Vi_num, Vo_num)*M_sec,'b','LineWidth',2)
plot(tf(phi_num,fs_num), iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(3), phi_num, fs_num, Vi_num, Vo_num)*M_sec,'g')
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