clear; close all; clc;

addpath('utils')
addpath('utils_transf')


%%

% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

color1 = [0.045 0.245 0.745]; % blue
color2 = [0.635 0.635 0.635]; % gray

syms fs Vi d Po dt t real positive
syms a1 b1 z1 a2 b2 z2 real
syms xa1 xb1 xz1 xa2 xb2 xz2 real
syms phi real
syms Ld1 Ld2 n Lm real positive
syms L1 L2 Ldab M real positive
syms fs Vi d dt real positive
syms phi real

Vo = Vi*d;
Ts = 1/fs;

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
n_num = 5/9;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = [700e-6, 1.4e-3 10e-3];
phi_num = [deg2rad(60)];  %[MUDAR]
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);


syms L1 L2 Ldab M real positive
syms fs Vi d dt real positive
syms phi real
syms Ld1 Ld2 n Lm real positive
    
Vo = Vi*d;
Ts = 1/fs;

x = sym('x_', [6,1], 'real');
dx = sym('dx_', [6,1], 'real');
eq = sym('eq_', [4,1], 'real');
u = sym('u_', [6,1], 'real');
s = sym('s_', [6,1], 'real');
    
% entradas
u(1:3) = s(1:3)*Vi;
u(4:6) = s(4:6)*Vo;
    
%corrente no primario do trafo
il = x(1:2);
il = [il; -sum(il)];
dil = dx(1:2);
dil = [dil; -sum(dil)];
    
%corrente no secundario do trafo
iL = x(4:5);
iL = [iL; -sum(iL)];
diL = dx(4:5);
diL = [diL; -sum(diL)];

%corrente no Ldab
ild = il - [il(3); il(1:2)];
dild = dil - [dil(3); dil(1:2)];

%definicao das tensoes do indutor acoplado
vP = + L1*dil + M*diL;
vS = + M*dil + L2*diL;

%definicao das tensoes nos indutores do dab
vLdab = Ldab*dild;

Td = [1 -1 0; 0 1 -1; -1 0 1];

%malhas
m_p = Td*u(1:3) == Td*vLdab + vP; %malha primario
m_s = Td*u(4:6) == Td*vS; %malha secundario

%usa 4 malhas
eq(1:2) = m_p(1:2);
eq(3:4) = m_s(1:2);

%% completar ilc e iLC
dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';

dxs(3) = -sum(dxs(1) + dxs(2));
dxs(6) = -sum(dxs(4) + dxs(5));
dxs = simplify(dxs);
Ms = equationsToMatrix(dxs, [x;s]);

As = Ms(:,1:6);
Bs = Ms(:,7:12);

%%


%% converter alpha beta
Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
Tclm = Tcl(1:2,:); %ignorar nivel zero
Tclxz = kron(eye(2), Tcl);
Tclx = kron(eye(2), Tclm);

Axx = simplify(Tclxz*As*inv(Tclxz));
Bxx = simplify(Tclxz*Bs*inv(Tclxz)); %a pinv deixa uns numeros feios
Ax = [Axx(:,1:2),Axx(:,4:5)];
Bx = [Bxx(:,1:2),Bxx(:,4:5)];
A = [Ax(1:2,:);Ax(4:5,:)];
B = [Bx(1:2,:);Bx(4:5,:)];
    
%% equivalent circuit


% syms dxa1 dxb1 dxa2 dxb2 real 
% syms xa1 xb1 xa2 xb2 real 
% syms ua1 ub1 ua2 ub2 real 
% 
% dxc = [dxa1 dxb1 n*dxa2 n*dxb2]';
% xc = [xa1 xb1 n*xa2 n*xb2]';
% uc = [ua1 ub1 ua2/n ub2/n]';
% 
% eq_cl = dxc == A*xc+B*uc;
% eq_cl = simplify(eq_cl);
% 
% eq_dx = struct2array(solve(eq_cl, [dxa1 dxb1 dxa2 dxb2],'Real',true,'IgnoreAnalyticConstraints',true))';
% eq_u = struct2array(solve(eq_cl, [ua1 ub1 ua2 ub2],'Real',true,'IgnoreAnalyticConstraints',true))';
%     
% M_dx = equationsToMatrix(eq_dx, [ua1 ub1 ua2 ub2]);
% 
% eigg = eig(M_dx)
% eigg = simplify(eigg)
% 
% 
% struct2array(solve(eq_cl, [dxa1 dxb1 dxa2 dxb2],'Real',true,'IgnoreAnalyticConstraints',true))'
% 
% 
% [V, J] = jordan(B)
% 
% J = simplify(J);
% 
% fJ = matlabFunction(J, 'vars', {L1,L2,Ldab,M});
% Jr = real(J)
% Ji = imag(J)
% Jn(1,:) = J(1,:)+J(3,:);
% Jn(2,:) = J(2,:)+J(4,:);
% Jn(3,:) = J(1,:)-J(3,:);
% Jn(4,:) = J(2,:)-J(4,:);
% Jn = simplify(Jn)
% 
% fJ = matlabFunction(J, 'vars', {L1,L2,Ldab,M});
% 
% fJ(L1_num(1), L2_num(1), Ldab_num, M_num(1))

%%
% T = sym('T', [2,2], 'real');
% syms a b c d real
% clear eqx
% lamb = lr + 1i*li;
% expand(simplify(a*d*lamb - b*c*lamb'))
% eqx(4) = det(T) == 1;

% solve(eqx, T)




% 
% 
% v1 = [1 1*1i 0];
% v2 = [1 0 1];
% v3 = [0 1 1];
% Mmm = [v1;v2;v3]
% rank(Mmm)

%% obter funcao de comutacao, depende de phi

[sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi);
sf = sf+0.5;
scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2

ts = simplify(ang)*Ts/(2*pi);
tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos

%% Obter valores de regime permanente

g = simplify(expm(A*dt));
h = B*dt;
xcl = sym('x_', [4,length(ts)], 'real');
x0 = sym('x0_',[4,1], 'real');

xcl(:,1) = x0;
for i=1:length(ts)-1
    dts = (ts(i+1)-ts(i));
    xcl(:,i+1) = simplify(subs(g, dt, dts)*xcl(:,i) + subs(h, dt, dts)*scl(:,i));
end

x0x = struct2array(solve(xcl(:,1) == -xcl(:,7), x0)).';
x0s = simplify(pinv(Tclx)*subs(xcl, x0, x0x));
dx0s = pinv(Tclx)*B*scl; %derivadas dos estados, indutor e trafo secundario

%% corrente rms e mean dos estados
[x_rms,~] = rms_and_mean(dx0s,x0s,ts,(1:12),(1:12));
IL_rms = x_rms(1);
Itrfsec_rms = x_rms(4);

f_IL_rms = matlabFunction(IL_rms, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
f_Itrfsec_rms = matlabFunction(Itrfsec_rms, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});

f_IL_rms(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)
f_Itrfsec_rms(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)

%% 
%obtido os valores de corrente na bobina do primario e secundario
%eh obtido todas as outras corrente

%corrente no Ldab
for ii=1:length(x0s)-1
    idab(:,ii) = subs(ild(1),x(1:3),x0s(1:3,ii));
    didab(:,ii) = subs(dild(1),dx(1:3),dx0s(1:3,ii));
end

hb = idab;
d_hb = didab;
HB = -x0s(4:6,:);
d_HB = -dx0s(4:6,:);

%% Corrente de entrada
target = [1;0;0];
[etapas] = pega_etapa(sf_p,target);
[irm,ime] = rms_and_mean(d_hb(1,:),hb(1,:),ts,etapas,etapas);
% f_ime = matlabFunction(ime, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
% f_irm = matlabFunction(irm, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
% 
% f_ime(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)
% f_irm(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)

%% Corrente de saida
target = [1;0;0];
[etapas] = pega_etapa(sf_s,target);
[iRM,iME] = rms_and_mean(d_HB(1,:),HB(1,:),ts,etapas,etapas);
% f_iME = matlabFunction(iME, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
% f_iRM = matlabFunction(iRM, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
% f_iME(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)
% f_iRM(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)

%% corrente de comutacao
Ip = hb(1,1);
Is = -HB(1,sec_switch);

%% funcao das corrente
fx0s = matlabFunction(x0s, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
fidab = matlabFunction(idab, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
fhb = matlabFunction(hb, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});
fHB = matlabFunction(HB, 'vars', {Ldab,L1,L2,M,phi,fs,Vi,d});

%%
ilL_num = fx0s(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num);
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
hold on
plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
plot(tf(phi_num,fs_num), ilL_num(1,:)','Color',jetcustom(1,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), ilL_num(4,:)','Color',jetcustom(2,:),'LineWidth',1.5)
hold off
grid on
grid minor
legend({'$i_{a}$[A]','$i_{A}$[A]'},'Location','best','FontSize', 14)
xlim([0 1/fs_num])
hold off
grid on
grid minor
set(gca, 'FontSize', 20)

%% plotar corrente

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

idab_num = fidab(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num);
idab_num = [idab_num, idab_num(:,1)];
figure
hold on
plot(tf(phi_num,fs_num), idab_num(1,:)','Color',jetcustom(2,:),'LineWidth',1.5)
hold off
xlim([0 1/fs_num])
hold off
grid on
grid minor
set(gca, 'FontSize', 20)





%% plotar corrente nos meia ponte

hb_num = fhb(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num);
hb_num = [hb_num, hb_num(:,1)];
HB_num = fHB(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num);


cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
hold on
plot(tf(phi_num,fs_num), hb_num(1,:)','Color',jetcustom(1,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), HB_num(1,:)','Color',jetcustom(2,:),'LineWidth',1.5)
hold off
xlim([0 1/fs_num])
hold off
grid on
grid minor
legend({'$i_{hb}$[A]','$i_{HB}$[A]'},'Location','best','FontSize', 14)
set(gca, 'FontSize', 20)


%%

fp = matlabFunction(ime*Vi,'vars',{Ldab,L1,L2,M,phi,fs,Vi,d});

figure
hold on
fplot(phi*180/pi,fp(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi, fs_num, Vi_num, d_num),[deg2rad(-30) 0],'Color',jetcustom(1,:),'LineWidth',1.5)

hold off
