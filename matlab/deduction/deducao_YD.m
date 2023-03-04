clear; close all; clc;

syms L1 L2 Ldab Lm real
syms n real

x = sym('x_', [6,1], 'real');
dx = sym('dx_', [6,1], 'real');
eq = sym('eq_', [6,1], 'real');
u = sym('u_', [6,1], 'real');

iL = x(1:2); %iLa e iLb
iL = [iL; -sum(iL)];
diL = dx(1:2);
diL = [diL; -sum(diL)];

iS = x(4:5); %iLA e iLB
iS = [iS; -sum(iS)];
diS = dx(4:5);
diS = [diS; -sum(diS)];

iD = iL;
diD = diL;

iM = iL - iS*n;
diM = diL - diS*n;

%definicao das tensoes
vM = Lm*diM;
vLdab = Ldab*diL;
vL1 = L1*diD;
vL2 = L2*diS;



Tdiff = [1 -1 0; 0 1 -1; -1 0 1];


eq(1:3) = Tdiff*(u(1:3) - vLdab - vM/n - vL1)  == 0;
eq(4:6) = Tdiff*u(4:6) == [vL2(2:3);vL2(1)] + [vM(2:3);vM(1)];


dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';

dxs(3) = -dxs(1)-dxs(2);
dxs(6) = -dxs(4)-dxs(5);

% Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
% Tclm = Tcl(1:2,:); %ignorar nivel zero
% Tclx = kron(eye(2), Tcl);

% dxx = Tclx.*dxs;
%%
M = equationsToMatrix(dxs, [x; s]);
A = M(:,1:2);
B = M(:,3:4); %derivadas de iLdab e iLd
    
%     A = kron(A, eye(2));
%     B = kron(B, eye(2));







%%


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
derivadas = pinv(Tclx)*B*scl; %derivadas dos estados, indutor e trafo secundario


%%

% solve(dxss, dx,'Real',true,'IgnoreAnalyticConstraints',true)
