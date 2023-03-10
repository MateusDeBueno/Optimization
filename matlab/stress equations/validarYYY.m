clear; close all; clc;

% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

addpath('utils');
addpath('utils_transf');
addpath('data');

syms phi real
syms fs Vi d Ldab Ld1 Ld2 Lm n Po dt t real positive

Vo = d*Vi;

% Parametros para validacao numerica
Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 10e-6;
Ld2_num = 6e-6;
Lm_num = [700e-6, 1.4e-3 10e-3];
phi_num = deg2rad(50);  %[MUDAR]
n_num = 5/9;
Po_num = 2000;
pi_num = 3.141592653589793;

variaveis =         [phi, fs, Vi, d, Ldab, Ld1, Ld2, Lm, n, Po];
ponto_de_operacao = [phi_num, fs_num, Vi_num, d_num, Ldab_num, Ld1_num, Ld2_num, Lm_num(1), n_num, Po_num];

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


solve(eq, s)

M = equationsToMatrix(dxs, [x; s]);
A = M(:,1:2);
B = M(:,3:4); %derivadas de iLdab e iLd
A = kron(A, eye(2));
B = kron(B, eye(2));

% syms Mc L1 L2 real positive
% old = [Ld1, Ld2, Ld2 + n*n*Lm];
% new = [Mc L1 L2];

%% Obter valores de regime permanente
g = simplify(expm(A*dt));
h = B*dt;
x = sym('x_', [4,length(ts)], 'real');
x0 = sym('x0_',[4,1], 'real');

x(:,1) = x0;
for i=1:length(ts)-1
    dts = (ts(i+1)-ts(i));
    x(:,i+1) = simplify(subs(g, dt, dts)*x(:,i) + subs(h, dt, dts)*scl(:,i));
%     x(:,i+1) = simplify(subs(h, dt, dts)*scl(:,i));
end

x0x = struct2array(solve(x(:,1) == -x(:,7), x0)).';
x0s = simplify(pinv(Tclx)*subs(x, x0, x0x));
iSw = matlabFunction(x0s.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, d});
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

%%  Validação dos Pulsos de Comutação
color1 = [0.045 0.245 0.745]; % blue
color2 = [0.635 0.635 0.635]; % gray

figure
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on
for ii=1:6
    if (mod(ii, 2) == 1)
        stairs(tf(phi_num, fs_num)*1e6, ([sf(ii,:) sf(ii,1)]'+2*ii),'-','Color',jetcustom(ceil(ii/2),:),'LineWidth',1.5)
    else 
        stairs(tf(phi_num, fs_num)*1e6, ([sf(ii,:) sf(ii,1)]'+2*ii),'--','Color',jetcustom(ceil(ii/2),:),'LineWidth',1.5)
    end
end
hold off
grid on
grid minor
xlim([0 1/fs_num]*1e6)
ylabel('$s_{abc,ABC}$')
xlabel('t [$\mu s$]')
set(gca, 'FontSize', 20)
legend({'$s_a$','$s_a$','$s_a$','$s_a$','$s_a$','$s_a$'},'Location','northwest','FontSize', 14,'NumColumns',3)



%%
% figure(1)
% subplot(2,1,1)
% stairs(tf(phi_num, fs_num) * 1e6, ([sf sf(:,1)]'+[0 -2 -4 -6 -8 -10]), 'LineWidth', 2)
% xlim([0 1/fs_num]*1e6)
% ylabel('$s_{abc,ABC}$')
% xlabel('t [$\mus$]')
% legend('Primary', 'Secondary')
% grid on
% subplot(2,1,2)
% stairs(tf(phi_num, fs_num) * 1e6, ([scl scl(:,1)]'+[0 -2 -4 -6]), 'LineWidth', 2)
% xlim([0 1/fs_num]*1e6)
% ylabel('s_{\alpha,\beta}')
% xlabel('t [\mus]')
% legend('Primary', 'Secondary')
% grid on



%% plots
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
M_sec = [[1 0 0];zeros(5,3)];
is1 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, d_num)*M_sec;
is2 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(2), phi_num, fs_num, Vi_num, d_num)*M_sec;
is3 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(3), phi_num, fs_num, Vi_num, d_num)*M_sec;
figure
hold on
plot(tf(phi_num,fs_num), is1(:,1)','Color',jetcustom(1,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), is2(:,1)','Color',jetcustom(2,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), is3(:,1)','Color',jetcustom(3,:),'LineWidth',1.5)
hold off
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlim([0 1/fs_num])
xlabel('$t$\,[s]')
ylabel('$i\,[A]$')
set(gca, 'FontSize', 20)


%% plot corrente secundario
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
M_sec = [zeros(3,3);[1 0 0];zeros(2,3)];
is1 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, d_num)*M_sec;
is2 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(2), phi_num, fs_num, Vi_num, d_num)*M_sec;
is3 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(3), phi_num, fs_num, Vi_num, d_num)*M_sec;
figure
hold on
plot(tf(phi_num,fs_num), is1(:,1)','Color',jetcustom(1,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), is2(:,1)','Color',jetcustom(2,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), is3(:,1)','Color',jetcustom(3,:),'LineWidth',1.5)
hold off
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlim([0 1/fs_num])
xlabel('$t$\,[s]')
ylabel('$i\,[A]$')
set(gca, 'FontSize', 20)

%% plot corrente da fase A, entrada e saida
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;
M1 = [1; 0; 0; 0; 0; 0];
M2 = [0; 0; 0; -1; 0; 0];
is1 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, d_num)*M1;
is2 = iSw(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, d_num)*M2;
figure
hold on
plot(tf(phi_num,fs_num), is1(:,1)','Color',jetcustom(1,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), is2(:,1)','Color',jetcustom(2,:),'LineWidth',1.5)
hold off
grid on
grid minor
legend({'$i_a$','$i_A$'},'Location','best','FontSize', 14)
xlim([0 1/fs_num])
xlabel('$t$\,[s]')
ylabel('$i\,[A]$')
set(gca, 'FontSize', 20)


%% ZVS d por phi
d_limit_p = solve(Ip == 0, d);
fd_limit_p = matlabFunction(d_limit_p, 'vars', {Ld2, Lm, phi, n});

d_limit_s = solve(Is == 0, d);
fd_limit_s = matlabFunction(d_limit_s, 'vars', {Ldab, Ld1, Lm, phi, n});

figure
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on
fplot(phi*180/pi,fd_limit_p(Ld2_num, Lm_num(1), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(1,:),'LineWidth',1.5)
fplot(phi*180/pi,fd_limit_p(Ld2_num, Lm_num(2), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(2,:),'LineWidth',1.5)
fplot(phi*180/pi,fd_limit_p(Ld2_num, Lm_num(3), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(3,:),'LineWidth',1.5)

fplot(phi*180/pi,fd_limit_s(Ldab_num, Ld1_num, Lm_num(1), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(1,:),'LineWidth',1.5)
fplot(phi*180/pi,fd_limit_s(Ldab_num, Ld1_num, Lm_num(2), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(2,:),'LineWidth',1.5)
fplot(phi*180/pi,fd_limit_s(Ldab_num, Ld1_num, Lm_num(3), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(3,:),'LineWidth',1.5)
hold off
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlabel('$\phi\,[^{\circ}]$')
ylabel('$d$')
set(gca, 'FontSize', 20)

%%
% figure
% fplot(phi,fd_limit_p(Ld2_num, Lm_num(1), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(1,:),'LineWidth',1.5)
% fplot(phi,fd_limit_s(Ldab_num, Ld1_num, Lm_num(1), phi, n_num),[0 deg2rad(60)],'Color',jetcustom(1,:),'LineWidth',1.5)



figure
hold on
f1 = fd_limit_p(Ld2_num, Lm_num(1), phi, n_num);
f2 = fd_limit_s(Ldab_num, Ld1_num, Lm_num(1), phi, n_num);
fill_between_functions(f1, f2, jetcustom(1,:), [0 deg2rad(60)])
f1 = fd_limit_p(Ld2_num, Lm_num(2), phi, n_num);
f2 = fd_limit_s(Ldab_num, Ld1_num, Lm_num(2), phi, n_num);
fill_between_functions(f1, f2, jetcustom(2,:), [0 deg2rad(60)])
f1 = fd_limit_p(Ld2_num, Lm_num(3), phi, n_num);
f2 = fd_limit_s(Ldab_num, Ld1_num, Lm_num(3), phi, n_num);
fill_between_functions(f1, f2, jetcustom(3,:), [0 deg2rad(60)])
hold off
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlabel('$\phi\,[rad]$')
ylabel('$d$')
set(gca, 'FontSize', 20)

%% ZVS
% phi_limit = solve(Ip == 0, phi);
% Vo_limit = Vo == Po/Iout_med;
% Vo_limit = subs(Vo_limit,phi,phi_limit);
% Vo_limit = isolate(Vo_limit,Vo);
% f_Vo_limit = matlabFunction(rhs(Vo_limit), 'vars', {Ldab, n, Ld1, Ld2, Lm, Po, fs, Vi});
% 
% figure
% fplot(Po,f_Vo_limit(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), Po, fs_num, Vi_num),[0 2000])
% 
% 

%%


%% calculo da variacao no Ldab
% 
% Vdab = Ldab*derivadas(1,:);
% fVdab = matlabFunction(Vdab.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
% fder = matlabFunction(derivadas(1,:).', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
% 
% time1 = tf(phi_num,fs_num);
% timex = kron(time1(1:end),ones(2));
% timey = timex(1,:);
% Vdab1 = fVdab(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, d_num);
% Vdab1 = [Vdab1; Vdab1(1)]';
% Vdabx = kron(Vdab1(1:end),ones(2));
% Vdaby = Vdabx(1,:);
% 
% figure
% plot(timey, circshift(Vdaby,[0,1]),'LineWidth',2)
% grid on
% grid minor
% xlim([0 1/fs_num])
% 
% 
% 
% %% calculo da variacao do fluxo no transformador (olhando pelo secundario)
% 
% Tm = kron([1 0;1 -1], eye(2)); %matriz de transformacao do circuito monofasico para fornecer
% %converter o segundo estado para corrente na magnetizando do primario
% 
% % Tm*B*scl
% 
% 
% Am = simplify(Tm*A*Tm^(-1));
% Bm = simplify(Tm*B);
% 
% x0m = pinv(Tclx)*Tclx*x0s;
% 
% % x0ma = x0s(1,:)-n*x0s(4,:);
% 
% x0m = x0s(1:3,:) - n*x0s(4:6, :)
% 
% iSm = matlabFunction(x0m.', 'vars', {Ldab, n, Ld1, Ld2, Lm, phi, fs, Vi, Vo});
% 
% iSw
% 
% x0s
% 
% 
% figure
% hold on
% plot(tf(phi_num,fs_num), iSm(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num(1), phi_num, fs_num, Vi_num, d_num)*eye(3,3),'r','LineWidth',2)
% hold off
% grid on
% grid minor
% xlim([0 1/fs_num])
% 
% derivadas_Vm = pinv(Tclx)*Bm*scl;
% derivadas_Vm(1,:).*ts(1:end-1)
% 
% 
% 
% 
% syms Vsec