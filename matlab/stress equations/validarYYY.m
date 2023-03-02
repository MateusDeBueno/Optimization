clear; close all; clc;
syms phi real
syms fs Vi Vo Ldab Ld1 Ld2 Lm n Po dt t real positive

% Parametros para validacao numerica
Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 10e-6;
Ld2_num = 5e-6;
Lm_num = [700e-6, 1e-3, 1000e-3];
phi_num = deg2rad(50);  %[MUDAR]
n_num = 5/9;
Vo_num = d_num*Vi_num;
Po_num = 2000;
pi_num = 3.141592653589793;

variaveis =         [phi, fs, Vi, Vo, Ldab, Ld1, Ld2, Lm, n, Po];
ponto_de_operacao = [phi_num, fs_num, Vi_num, Vo_num, Ldab_num, Ld1_num, Ld2_num, Lm_num(1), n_num, Po_num];



Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
Tclm = Tcl(1:2,:); %ignorar nivel zero
Tclx = kron(eye(2), Tclm);
    
Ts = 1/fs;
    
[sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi);

scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
ts = simplify(ang)*Ts/(2*pi);
tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos

%% Resolvendo circuito monofasico equivalente [MUDAR]
L = [Ldab+Ld1; Ld2/n^2];
Udc = [Vi; Vo];
u = Udc.*[1; 1/n];
    
x = sym('x_', [2,1], 'real');
dx = sym('dx_', [2,1], 'real');
eq = sym('eq', [2,1], 'real');
s = sym('s', [2,1], 'real');

Vm = Lm*(dx(1)-dx(2)*n);
VL1 = L(1)*dx(1);
VL2 = L(2)*dx(2)*n;

eq(1) = u(1)*s(1) == VL1 + Vm;
eq(2) = u(2)*s(2) == Vm - VL2;

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
[x_rms,~] = rms_and_mean(derivadas,x0s,ts,(1:12),(1:12));
IL_rms = x_rms(1);
Itrfsec_rms = x_rms(4);

%% Fourier dos estados
% IL_rms_c_k = ck_fourier(ts,derivadas(1,:),x0s(1,:)); %DESCOMENTAR
% Itrf_sec_rms_c_k = ck_fourier(ts,derivadas(4,:),x0s(4,:)); %DESCOMENTAR
    
%% Corrente do half bridge do primario [MUDAR]
derivadas_hb_p = derivadas(1:3,:);
pts_inics_hb_p = x0s(1:3,:);
    
%% Corrente do half bridge do secundario [MUDAR]
derivadas_hb_s = derivadas(4:6,:);
pts_inics_hb_s = x0s(4:6,:);

%% Corrente de entrada
target = [1;0;0];
[etapas] = pega_etapa(sf_p,target);
[Iin_rms,Iin_med] = rms_and_mean(derivadas_hb_p(1,:),pts_inics_hb_p(1,:),ts,etapas,etapas);
    
%% Corrente de saida
target = [1;0;0];
[etapas] = pega_etapa(sf_s,target);
[Iout_rms,Iout_med] = rms_and_mean(derivadas_hb_s(1,:),pts_inics_hb_s(1,:),ts,etapas,etapas);
    
power_equation = Iout_med==Po/Vo;
[phi_sol,~,phi_cond] = solve(power_equation,phi,'ReturnConditions',true);

%% Corrente nas chaves
[Isw_p_rms,~] = rms_and_mean(derivadas_hb_p(1,:),pts_inics_hb_p(1,:),ts,(1:6),(1:12));

[Isw_s_rms,~] = rms_and_mean(derivadas_hb_s(1,:),pts_inics_hb_s(1,:),ts,(1:6),(1:12));

%% corrente de comutacao
%primario
Ip = pts_inics_hb_p(1,1);

%secundario
Is = -pts_inics_hb_s(1,sec_switch);

%% output
output = [IL_rms,Itrfsec_rms,Iin_med,Iin_rms,Iout_med,Iout_rms,Isw_p_rms,Isw_s_rms,Ip,Is];
% output_harmonic = [IL_rms_c_k,Itrf_sec_rms_c_k];
output_phi = [phi_sol;phi_cond];

















%%  Validação dos Pulsos de Comutação
figure(1)
subplot(2,1,1)
stairs(tf(phi_num, fs_num) * 1e6, ([sf sf(:,1)]'+[0 -2 -4 -6 -8 -10]), 'LineWidth', 2)
xlim([0 1/fs_num]*1e6)
ylabel('s_{abc,ABC}')
xlabel('t [\mus]')
legend('Primary', 'Secondary')
grid on
subplot(2,1,2)
stairs(tf(phi_num, fs_num) * 1e6, ([scl scl(:,1)]'+[0 -2 -4 -6]), 'LineWidth', 2)
xlim([0 1/fs_num]*1e6)
ylabel('s_{\alpha,\beta}')
xlabel('t [\mus]')
legend('Primary', 'Secondary')
grid on


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



%% calculo da variacao no Ldab

Vdab = Ldab*derivadas(1,:);
fVdab = matlabFunction(Vdab.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
fder = matlabFunction(derivadas(1,:).', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

time1 = tf(phi_num,fs_num);
timex = kron(time1(1:end),ones(2));
timey = timex(1,:);
Vdab1 = fVdab(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num);
Vdab1 = [Vdab1; Vdab1(1)]';
Vdabx = kron(Vdab1(1:end),ones(2));
Vdaby = Vdabx(1,:);

figure
plot(timey, circshift(Vdaby,[0,1]),'LineWidth',2)
grid on
grid minor
xlim([0 1/fs_num])



%% calculo da variacao do fluxo no transformador (olhando pelo secundario)

Tm = kron([1 0;1 -1], eye(2)); %matriz de transformacao do circuito monofasico para fornecer
%converter o segundo estado para corrente na magnetizando do primario

% Tm*B*scl


Am = simplify(Tm*A*Tm^(-1));
Bm = simplify(Tm*B);

x0m = pinv(Tclx)*Tclx*x0s;

% x0ma = x0s(1,:)-n*x0s(4,:);

x0m = x0s(1:3,:) - n*x0s(4:6, :)

iSm = matlabFunction(x0m.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});

iSw

x0s


figure
hold on
plot(tf(phi_num,fs_num), iSm(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, Vo_num)*eye(3,3),'r','LineWidth',2)
hold off
grid on
grid minor
xlim([0 1/fs_num])

derivadas_Vm = pinv(Tclx)*Bm*scl;
derivadas_Vm(1,:).*ts(1:end-1)




syms Vsec