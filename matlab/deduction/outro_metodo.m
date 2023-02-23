clear; close all; clc;

syms phi fs Vi Vo real

Tcl = 1/3*[2 -1 -1; 0 sqrt(3) -sqrt(3)];
Tcli = pinv(Tcl);

%estados de comutacao de um inversor 3phase, ponto medio no meio do barramento
s = -0.5+[1 0 1; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]';
% s = 0.5-[1 0 1; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]';

s = kron(s, [1, 1]); %abc
s = [s; s(:, end) s(:,1:end-1)]; %a1b1c1 to a2b2c2
s = [s s(:,1)]; %poem ultima na ultima coluna a primeira

T = 1/fs;

%tempo de comutacao, valido para phi < 60 graus
t = [0, phi, pi/3, pi/3+phi, 2*pi/3, 2*pi/3+phi, pi, pi+phi, pi+pi/3, pi+pi/3+phi, pi+2*pi/3, pi+2*pi/3+phi, 2*pi]*T/(2*pi); 
tf = matlabFunction(t, 'vars', {phi, fs});

%transformada de clarke dos 4 estados de comutacao
Tclx = kron(Tcl, eye(2)); 
Tclxi = pinv(Tclx); %e inversa

%estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
scl = kron(eye(2), Tcl)*s;


%%
syms Ldab Ld Lm n real

x = sym('x_', [2,1], 'real');
dx = sym('dx_', [2,1], 'real');
equation = sym('eq', [2,1], 'real');
s = sym('s', [2,1], 'real');

equation(1) = s(1)*Vi == (Ldab + Ld)*dx(1) + Lm*(dx(1)+dx(2)/n);
equation(2) = Lm*(dx(1)+dx(2)/n)*n == Ld*dx(2) + s(2)*Vo;
dxs = struct2array(solve(equation, dx)).';
M = equationsToMatrix(dxs, [x; s]);
A = M(:,1:2);
B = M(:,3:4); %devidas de iLdab e iLd

% A = kron(A, eye(2));
B = kron(B, eye(2));

% B = blkdiag(B, B); %devidas de iLdab e iLd para alpha e beta

x0 = sym('x0_', [4,1], 'real');
x = sym('x_', [4,length(t)], 'real');

x(:,1) = x0;
vs = [1 2 3 4 5 6 7 8 9 10 11 12];
% vs = [1 3 2 4 5 4 8 7 10 9 12 11];
for i=1:length(t)-1
    dt = t(i+1) - t(i); %aqui tem uma simplificacao para somente retas
    x(:,i+1) = simplify(x(:,i) + dt*B*scl(:,vs(i))); %tempo x derivada x estados
end

%%
eq = x(:,1) == -x(:,7);
x0x = struct2array( solve(eq, x0)).';
xs = subs(x, x0, x0x);
xt = Tclxi*xs;


xsf = matlabFunction(xs, 'vars', {Ldab,Ld,Lm,Vi,Vo,fs,n,phi});
xtf = matlabFunction(xt, 'vars', {Ldab,Ld,Lm,Vi,Vo,fs,n,phi});


%%

tx = tf(deg2rad(50), 100e3);
xsw = xsf(60e-6, 2e-6, 0.7e-3, 400, 400, 100e3, 5/9, deg2rad(50));
xtw = xtf(60e-6, 2e-6, 0.7e-3, 400, 400, 100e3, 5/9, deg2rad(50));


figure
plot(tx, xsw(1:4,:)')

figure
plot(tx, xtw(1:3,:)')

