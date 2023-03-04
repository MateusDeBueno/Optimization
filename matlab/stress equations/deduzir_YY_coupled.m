clear; close all; clc;

addpath('utils')
addpath('utils_transf')
syms fs Vi d Po dt t real positive
syms a1 b1 z1 a2 b2 z2 real
syms xa1 xb1 xz1 xa2 xb2 xz2 real
syms phi real
syms Ld1 Ld2 n Lm real positive

Vo = Vi*d;
Ts = 1/fs;

% Parametros para validacao numerica
Vi_num = 400;
d_num = 0.5;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
Ld2_num = 0.6e-6;
Lm_num = 700e-6;
phi_num = deg2rad(50);  %[MUDAR]
n_num = 5/9;
Po_num = 2000;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num*L2_num);

syms L1 L2 Ldab M real positive

x = sym('x_', [6,1], 'real');
dx = sym('dx_', [6,1], 'real');
eq = sym('eq_', [4,1], 'real');
u = sym('u_', [6,1], 'real');
s = sym('s_', [6,1], 'real');

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
%definicao das tensoes do indutor acoplado
vP = + L1*dil - M*diL;
vS = - M*dil + L2*diL;
%definicao das tensoes nos indutores do dab
vLdab = Ldab*dil;


Tdiff = [1 -1 0; 0 1 -1; -1 0 1];

m_p = -Tdiff*u(1:3) + Tdiff*vLdab + Tdiff*vP  == 0; %malha primario
m_s = -Tdiff*u(4:6) - Tdiff*vS == 0; %malha secundario



eq(1:2) = m_p(1:2);
eq(3:4) = m_s(1:2);

%%
dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';

dxs(3) = -sum(dxs(1) + dxs(2));
dxs(6) = -sum(dxs(4) + dxs(5));
dxs = simplify(dxs);
Ms = equationsToMatrix(dxs, [x;s]);

As = Ms(:,1:6);
Bs = Ms(:,7:12);

%% converter alpha beta
Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
Tclm = Tcl(1:2,:); %ignorar nivel zero
Tclxz = kron(eye(2), Tcl);
Tclx = kron(eye(2), Tclm);

% A = simplify(Tclx*As*pinv(Tclx));
% B = simplify(Tclx*Bs*pinv(Tclx));
Axx = simplify(Tclxz*As*inv(Tclxz));
Bxx = simplify(Tclxz*Bs*inv(Tclxz)); %a pinv deixa uns numeros feios
Ax = [Axx(:,1:2),Axx(:,4:5)];
Bx = [Bxx(:,1:2),Bxx(:,4:5)];
A = [Ax(1:2,:);Ax(4:5,:)];
B = [Bx(1:2,:);Bx(4:5,:)];
%%

[sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi);
sf = sf+0.5;
scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2

ts = simplify(ang)*Ts/(2*pi);
tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos

% syms real
% ucl = sym('ucl_', [4,1], 'real');
% dcl = sym('dcl_', [4,1], 'real');
% 
% ducl = B*ucl == dcl;
% 
% solve(ducl, dcl)
% 
% 
% solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)
% dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).'

%%

figure
subplot(2,1,1)
stairs(tf(phi_num, fs_num) * 1e6, ([sf sf(:,1)]'+[0 -2 -4 -6 -8 -10]), 'LineWidth', 2)
xlim([0 1/fs_num]*1e6)
grid on
subplot(2,1,2)
stairs(tf(phi_num, fs_num) * 1e6, ([scl scl(:,1)]'+[0 -2 -4 -6]), 'LineWidth', 2)
xlim([0 1/fs_num]*1e6)
grid on

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
iSw = matlabFunction(x0s.', 'vars', {Ldab, L1, L2, M, phi, fs, Vi, d});

derivadas = pinv(Tclx)*B*scl; %derivadas dos estados, indutor e trafo secundario




%%
color1 = [0.045 0.245 0.745]; % blue
color2 = [0.635 0.635 0.635]; % gray

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;
M1 = [1; 0; 0; 0; 0; 0];
M2 = [0; 0; 0; 1; 0; 0];
is1 = iSw(Ldab_num, L1_num, L2_num, M_num, phi_num, fs_num, Vi_num, d_num)*M1;
is2 = iSw(Ldab_num, L1_num, L2_num, M_num, phi_num, fs_num, Vi_num, d_num)*M2;
figure
hold on
plot(tf(phi_num,fs_num), is1(:,1)','Color',jetcustom(1,:),'LineWidth',1.5)
plot(tf(phi_num,fs_num), is2(:,1)','Color',jetcustom(2,:),'LineWidth',1.5)
hold off
grid on
grid minor
legend({'i_a','i_A'},'Location','best','FontSize', 14)
xlim([0 1/fs_num])
hold off
grid on
grid minor
xlim([0 1/fs_num])