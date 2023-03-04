clear; close all; clc;

addpath('utils')
addpath('utils_transf')

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
pi_num = 3.141592653589793;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num*L2_num);




syms L1 L2 Ldab M real positive

x = sym('x_', [6,1], 'real');
dx = sym('dx_', [6,1], 'real');
eq = sym('eq_', [4,1], 'real');
u = sym('u_', [6,1], 'real');

iLp = x(1:2);
iLp = [iLp; -sum(iLp)];
diLp = dx(1:2);
diLp = [diLp; -sum(diLp)];

iLs = x(4:5);
iLs = [iLs; -sum(iLs)];
diLs = dx(4:5);
diLs = [diLs; -sum(diLs)];

%corrente na bobina do primario
% ii = (iLp(1:2) - iLp(2:3))/3;
% ii = [ii; -sum(ii)];
% dii = (diLp(1:2) - diLp(2:3))/3;
% dii = [dii; -sum(dii)];


%definicao das tensoes do indutor acoplado
vP = + L1*diLp - M*diLs;
vS = - M*diLp + L2*diLs;

%definicao das tensoes nos indutores
vLdab = Ldab*diLp;


Tdiff = [1 -1 0; 0 1 -1; -1 0 1];
m_p = -Tdiff*u(1:3) + vLdab + vP  == 0;
m_s = -Tdiff*u(4:6) - Tdiff*vS == 0;

eq(1:2) = m_p(1:2);
eq(3:4) = m_s(1:2);


dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';

% dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';

dxs(3) = -sum(dxs(1) + dxs(2));
dxs(6) = -sum(dxs(4) + dxs(5));
dxs = simplify(dxs);
Ms = equationsToMatrix(dxs, [x;u]);

As = Ms(:,1:6);
Bs = Ms(:,7:12);

%% converter alpha beta
Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
Tclx = kron(eye(2), Tcl);

Ax = Tclx*As*Tclx^-1;
Bx = Tclx*Bs*Tclx^-1;


syms a1 b1 z1 a2 b2 z2 real
syms xa1 xb1 xz1 xa2 xb2 xz2 real


ucl = [a1 b1 z1 a2 b2 z2]';
xcl = [xa1 xb1 xz1 xa2 xb2 xz2]';

dabz = simplify(Ax*xcl + Bx*ucl);