clear; close all; clc;

syms L1 L2 Ldab Lm real
syms n real

x = sym('x_', [4,1], 'real');
dx = sym('dx_', [4,1], 'real');
eq = sym('eq_', [4,1], 'real');
u = sym('u_', [6,1], 'real');

iL = x(1:2);
iL = [iL; -sum(iL)];
diL = dx(1:2);
diL = [diL; -sum(diL)];



iS = x(3:4);
iS = [iS; -sum(iS)];
diS = dx(3:4);
diS = [diS; -sum(diS)];

%%

iD = 





% iD = x(3:4);
% iD = [iD; -sum(iD)];
% diD = dx(3:4);
% diD = [diD; -sum(diD)];

iM = iD - n*iS;
diM = diD - n*diS;




%definicao das tensoes
vM = Lm*diM;
vLdab = Ldab*diL;
vL1 = L1*diD;
vL2 = L2*diS;

%%

Tdiff = [1 -1 0; 0 1 -1; -1 0 1];


eq(1:3) = Tdiff*u(1:3) - Tdiff*vLdab - vM/n - vL1  == 0
eq(4) = sum(vL1 + vM/n) == 0
eq(5:6) = Tdiff(1:2,:)*(u(4:6) - vL2 - vM) == 0


% eq(1) = -u(1) + vLdab(1) + vM(1)/n +vL1(1) - vLdab(2) + u(2) == 0;
% eq(2) = -u(2) + vLdab(2) + vM(2)/n +vL1(2) - vLdab(3) + u(3) == 0;
% eq(3) = -u(3) + vLdab(3) + vM(3)/n +vL1(3) - vLdab(1) + u(1) == 0;
% eq(5) = -u(4) + vL2(1) + vM(1) - vM(2) - vL2(2) + u(5) == 0;
% eq(6) = -u(5) + vL2(2) + vM(2) - vM(3) - vL2(3) + u(6) == 0;



dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';

dxss = subs(eq,dx(1:4),dxs(1:4));

% solve(dxss, dx,'Real',true,'IgnoreAnalyticConstraints',true)
