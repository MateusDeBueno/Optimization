clear; close all; clc;

color1 = [0.045 0.245 0.745]; % blue
color2 = [0.635 0.635 0.635]; % gray


% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

addpath('utils\')
addpath('utils_transf\')

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
n_num = 1;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = 700e-6;
phi_num = deg2rad(10);  %[MUDAR]
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

%% aqui comeca
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

%entradas
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
vP = + L1*dil + M*diL;
vS = + M*dil + L2*diL;

%matrix auxiliar
Td = [1 -1 0; 0 1 -1; -1 0 1];

%% definicoes primario %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%corrente no Ldab
ild = il - [il(3); il(1:2)];
dild = dil - [dil(3); dil(1:2)];
%corrente no HB
ihb = il - [il(3); il(1:2)];
dihb = dil - [dil(3); dil(1:2)];
%tensoes no Ldab
vLdab = Ldab*dild;
%malha primario
m_p = Td*u(1:3) == Td*vLdab + vP; 

%% definicoes secundario %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%corrente no HB
iHB = iL - [iL(3); iL(1:2)];
diHB = diL - [diL(3); diL(1:2)];
%malha secundario
m_s = Td*u(4:6) == vS; 

%% usa 4 malhas
eq(1:2) = m_p(1:2);
eq(3:4) = m_s(1:2);

%resolve
dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';

%% completar ilc e iLC
dxs(3) = -sum(dxs(1) + dxs(2));
dxs(6) = -sum(dxs(4) + dxs(5));
dxs = simplify(dxs);
Ms = equationsToMatrix(dxs, [x;s]);
As = Ms(:,1:6);
Bs = Ms(:,7:12);

%% converter alpha beta
Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
Tclm = Tcl(1:2,:); %ignorar nivel zero
Tclx = kron(eye(2), Tclm);
A = Tclx*As*pinv(Tclx);
B = Tclx*Bs*pinv(Tclx);

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
    
%% Corrente nos estados
[ilrm,~] = rms_and_mean(dx0s(1,:),x0s(1,:),ts,1:12,1:12);
[iLrm,~] = rms_and_mean(dx0s(4,:),x0s(4,:),ts,1:12,1:12);

%% Corrente hb
for ii=1:length(x0s)-1
    hb(:,ii) = subs(ihb(1),x,x0s(:,ii));
    dhb(:,ii) = subs(dihb(1),dx,dx0s(:,ii));
end

%% Corrente HB
for ii=1:length(x0s)-1
    HB(:,ii) = subs(iHB(1),x,x0s(:,ii));
    dHB(:,ii) = subs(diHB(1),dx,dx0s(:,ii));
end

%% Corrente Ldab
for ii=1:length(x0s)-1
    idab(:,ii) = subs(ild(1),x,x0s(:,ii));
    didab(:,ii) = subs(dild(1),dx,dx0s(:,ii));
end
[idrm,~] = rms_and_mean(didab(1,:),idab(1,:),ts,1:12,1:12);

%% Corrente nas chaves
[iSwPrm,~] = rms_and_mean(dhb(1,:),hb(1,:),ts,1:6,1:12);
[iSwSrm,~] = rms_and_mean(-dHB(1,:),-HB(1,:),ts,1:6,1:12);  

%% corrente de comutacao
Ip = hb(1,1);
Is = -HB(1,sec_switch);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pr = 200; %precision
vec_ph = min(intervalo):(max(intervalo)-min(intervalo))/pr:max(intervalo);

fx0s = matlabFunction(x0s, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fts = matlabFunction(ts, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

fidrm = matlabFunction(idrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
filrm = matlabFunction(ilrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiLrm = matlabFunction(iLrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwPrm = matlabFunction(iSwPrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwSrm = matlabFunction(iSwSrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

fIp = matlabFunction(Ip, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs = matlabFunction(Is, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

mag=1;

fIp(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1))
fIs(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1))

%%
equat = fx0s(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1));
equat_ts = fts(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1));

figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap);
jetcustom = cmap;
hold on
plot(equat_ts*1e6,equat(1,:),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(equat_ts*1e6,-equat(4,:), '-.','Color',jetcustom(2,:),'LineWidth',1.5)
hold off
xlim([0 max(equat_ts)*1e6])
grid on
grid minor
legend({'$i_l$','$-i_L$'},'Location','best','FontSize', 14)
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
