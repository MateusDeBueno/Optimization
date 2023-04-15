clear; close all; clc;
%%
color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue

addpath('dados_pratica')

% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end



addpath('utils')
addpath('utils_transf')
addpath('utils_loss')
addpath('data');

awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32

l.eq = load('YY.mat');
l.sw = load('f_fitted_off.mat');
l.sw = load('f_fitted_on.mat');
l.sw.Ronp = 90e-3;
l.sw.Rons = 90e-3;


l.L.a = 1.394;
l.L.b = 2.248;
l.L.kc = 2.448;
l.L.int_ki = integral(@(theta) abs(cos(theta)).^l.L.a,0,2*pi);
l.L.ki = l.L.kc/(2^(l.L.b-l.L.a)*(2*pi)^(l.L.a-1)*l.L.int_ki);
l.L.Ac = 421.3e-6;
l.L.Ve = 52.1e-6;
l.L.N = 14;
l.L.strand = 180;
l.L.awg = 38;
l.L.dl = awg_wires(l.L.awg)*1e-3; % [m]
l.L.MLT = 140*1e-3; % [m]
l.L.Nl = 1;

l.tr.a = 1.585;
l.tr.b = 2.756;
l.tr.kc = 0.402;
l.tr.int_ki = integral(@(theta) abs(cos(theta)).^l.tr.a,0,2*pi);
l.tr.ki = l.tr.kc/(2^(l.tr.b-l.tr.a)*(2*pi)^(l.tr.a-1)*l.tr.int_ki);
l.tr.Ac = 683e-6;
l.tr.Ve = 146.79e-6;
l.tr.Np = 9;
l.tr.Ns = 5;
l.tr.Sp = 2;
l.tr.Ss = 3;
l.tr.N = l.tr.Ns;
l.tr.strand = 180;
l.tr.awg = 38;
l.tr.dl = awg_wires(l.tr.awg)*1e-3; % [m]
l.tr.MLT = 230*1e-3; % [m]

l.sC.Rac100 = 4.3e-3;


l.pr.dt = 0e-9;
l.pr.phi = deg2rad(0);
l.pr.Vi = 400;
l.pr.d = 1;
l.pr.fs = 100e3;
l.pr.Ldab = 61.5e-6;
l.pr.Ld1 = 1.5e-6;
l.pr.n = l.tr.Ns/l.tr.Np;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 700e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);

Vi_num = l.pr.Vi;
d_num = l.pr.d;
fs_num = l.pr.fs;
Ldab_num = l.pr.Ldab;
Ld1_num = l.pr.Ld1;
n_num = l.pr.n;
Ld2_num = l.pr.Ld2;
Lm_num = l.pr.Lm;
phi_num = l.pr.phi;  %[MUDAR]
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);











%cria intervalo de angulo
intervalo = [0 deg2rad(60)];


%% aqui comeca

    syms L1 L2 Ldab M real positive
    syms fs Vi d dt real positive
    syms phi real
    
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
    ild = il;
    dild = dil;
    %corrente no HB
    ihb = il - [il(3); il(1:2)];
    dihb = dil - [dil(3); dil(1:2)];
    %tensoes no Ldab
    vLdab = Ldab*dild;
    %malha primario
    m_p = Td*u(1:3) == vLdab + vP; 
    
    %% definicoes secundario %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %corrente no HB
    iHB = -iL;
    diHB = -diL;
    %malha secundario
    m_s = Td*u(4:6) == Td*vS; 
    
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
    
    %fourier dos estados
    syms nn
%     ilrm_cn = ck_fourier(ts,dx0s(1,:),x0s(1,:)); %DESCOMENTAR
%     iLrm_cn = ck_fourier(ts,dx0s(4,:),x0s(4,:)); %DESCOMENTAR  
    
    %% Corrente hb
    for ii=1:length(x0s)-1
        hb(:,ii) = subs(ihb(1),x,x0s(:,ii));
        dhb(:,ii) = subs(dihb(1),dx,dx0s(:,ii));
    end
    [hbrm,~] = rms_and_mean(dhb(1,:),hb(1,:),ts,1:12,1:12);
    
    %% Corrente HB
    for ii=1:length(x0s)-1
        HB(:,ii) = subs(iHB(1),x,x0s(:,ii));
        dHB(:,ii) = subs(diHB(1),dx,dx0s(:,ii));
    end
    [HBrm,~] = rms_and_mean(dHB(1,:),HB(1,:),ts,1:12,1:12);
    
    %% Corrente de entrada
    target = [1;0;0];
    [etapas] = pega_etapa(sf_p,target);
    [iiRMS,iiME] = rms_and_mean(dhb(1,:),hb(1,:),ts,etapas,etapas);
    
    %% Corrente de saida
    target = [1;0;0];
    [etapas] = pega_etapa(sf_s,target);
    [ioRMS,ioME] = rms_and_mean(dHB(1,:),HB(1,:),ts,etapas,etapas);
    Pm = Vo*ioME;
    
    %% Corrente Ldab
    for ii=1:length(x0s)-1
        idab(:,ii) = subs(ild(1),x,x0s(:,ii));
        didab(:,ii) = subs(dild(1),dx,dx0s(:,ii));
    end
    [idrm,~] = rms_and_mean(didab(1,:),idab(1,:),ts,1:12,1:12);
%     idrm_cn = ck_fourier(ts,didab(1,:),idab(1,:)); %DESCOMENTAR
    
    %% Corrente nas chaves
    [iSwPrm,~] = rms_and_mean(dhb(1,:),hb(1,:),ts,1:6,1:12);
    [iSwSrm,~] = rms_and_mean(dHB(1,:),HB(1,:),ts,1:6,1:12);  
    
    %% corrente de comutacao
    Ip = hb(1,1);
    Is = -HB(1,sec_switch);
    
    %% inductor core loss
    syms ki_L b_L a_L Ac_L N_L Ve_L real positive
    
    wb_L = Ldab*idab;
    dwb_L = Ldab*didab;
    B_L = wb_L/(N_L*Ac_L);
    dB_L = dwb_L/(N_L*Ac_L);
    
    integrais = sym('integrais_',[1,length(ts)-1], 'real');
    for ii=1:length(ts)-1
        dts = (ts(ii+1)-ts(ii));
        integrais(ii) = dts*abs(dB_L(ii))^a_L;
    end
    integrais = simplify(integrais);
    sum_int = simplify(sum(integrais));
    
    Bpk_L = max(simplify(B_L));
    Bppk_L = 2*Bpk_L;
    
    Pv_L = simplify(ki_L*fs*sum_int)*Bppk_L^(b_L-a_L);
    P_core_L = Pv_L*Ve_L;
    
    %% secondary core loss
    syms ki_tr b_tr a_tr Ac_tr N_tr Ve_tr real positive
    
    for ii=1:length(x0s)-1
        vSS(ii) = subs(vS(1),dx,dx0s(:,ii));
    end
    vSS = simplify(vSS);
    
    dB_tr = vSS/(N_tr*Ac_tr);
    
    integrais = sym('integrais_',[1,length(ts)-1], 'real');
    B_tr = sym('B_tr_',[1,length(ts)-1], 'real');
    B_tr(1) = 0;
    for ii=1:length(ts)-1
        dts = (ts(ii+1)-ts(ii));
        integrais(ii) = dts*abs(dB_tr(ii))^a_tr;
        B_tr(ii+1) = B_tr(ii) + dB_tr(ii)*dts;
    end
    integrais = simplify(integrais);
    sum_int_tr = simplify(sum(integrais));
    B_tr = simplify(B_tr);
    
    Bppk_tr = max(B_tr)-min(B_tr);
    Bpk_tr = Bppk_tr/2;
    
    Pv_tr = simplify(ki_tr*fs*sum_int_tr)*Bppk_tr^(b_tr-a_tr);
    P_core_tr = Pv_tr*Ve_tr;
    
    
%%

fts = matlabFunction(ts, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fHB = matlabFunction(HB, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fhb = matlabFunction(hb, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fhbrm = matlabFunction(hbrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fHBrm = matlabFunction(HBrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fidab = matlabFunction(idab, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIp = matlabFunction(Ip, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs = matlabFunction(Is, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fPm = matlabFunction(Pm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fidrm = matlabFunction(idrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
filrm = matlabFunction(ilrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiLrm = matlabFunction(iLrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwPrm = matlabFunction(iSwPrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwSrm = matlabFunction(iSwSrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

fiiRMS = matlabFunction(iiRMS, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiiME = matlabFunction(iiME, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fioRMS = matlabFunction(ioRMS, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

fP_core_tr = matlabFunction(P_core_tr, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi,N_tr,Ac_tr,Ve_tr,a_tr,b_tr,ki_tr});
fBpk_tr = matlabFunction(Bpk_tr, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi,N_tr,Ac_tr});
fP_core_L = matlabFunction(P_core_L, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi,N_L,Ac_L,Ve_L,a_L,b_L,ki_L});
fBpk_L = matlabFunction(Bpk_L, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi,N_L,Ac_L});




fsum_int_tr = matlabFunction(sum_int_tr, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi,N_tr,Ac_tr,Ve_tr,a_tr,b_tr,ki_tr});



fhbrm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fHBrm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fIp(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fIs(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fPm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fidrm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
filrm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fiLrm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fiSwPrm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fiSwSrm(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)

fiiRMS(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fiiME(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)
fioRMS(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num)

fP_core_tr(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.tr.N,l.tr.Ac,l.tr.Ve,l.tr.a,l.tr.b,l.tr.ki)
fBpk_tr(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.tr.N,l.tr.Ac)
fP_core_L(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.L.N,l.L.Ac,l.L.Ve,l.L.a,l.L.b,l.L.ki)
fBpk_L(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.L.N,l.L.Ac)

Pse = l.tr.Ve*l.tr.kc*(l.pr.fs)^(l.tr.a)*fBpk_tr(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.tr.N,l.tr.Ac)^(l.tr.b);

%% criar uma porrada de funcoes
syms Po real

fidab = matlabFunction(idab, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fx0s = matlabFunction(x0s, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fts = matlabFunction(ts, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%equacao_limite
Ip_eq = Ip == 0;
Is_eq = Is == 0;
pot_eq = ioME*Vo == Po;

fpot_eq = matlabFunction(solve(pot_eq,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIp_eq = matlabFunction(solve(Ip_eq,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs_eq = matlabFunction(solve(Is_eq,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%correntes eficazes
fidrm = matlabFunction(idrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
filrm = matlabFunction(ilrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiLrm = matlabFunction(iLrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwPrm = matlabFunction(iSwPrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwSrm = matlabFunction(iSwSrm, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

fIp = matlabFunction(Ip, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs = matlabFunction(Is, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%% coisas para parte teorica
% corrente nos estados

trafo = 'DfY';


yy = fx0s(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num(1));
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num(1));

%corrente no primario do trafo
figure
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on

plot(xx*1e6,yy(1,:),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx*1e6,yy(2,:),'Color',jetcustom(2,:),'LineWidth',1.5)
plot(xx*1e6,yy(3,:),'Color',jetcustom(3,:),'LineWidth',1.5)

hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
legend({'$i_{a}$','$i_{b}$','$i_{c}$'},'Location','southeast','FontSize', 16)
file_name = append('figure\finalCap2\primary_current_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');


%corrente no secundario do trafo
figure
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on

plot(xx*1e6,yy(4,:),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx*1e6,yy(5,:),'Color',jetcustom(2,:),'LineWidth',1.5)
plot(xx*1e6,yy(6,:),'Color',jetcustom(3,:),'LineWidth',1.5)

hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
legend({'$i_{A}$','$i_{B}$','$i_{C}$'},'Location','southeast','FontSize', 16)
file_name = append('figure\finalCap2\secondary_current_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

% corrente no Ldab
yy_s = fx0s(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num(1));
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num(1));
il = yy_s(1:3,:);
yy = il - [il(3,:); il(1:2,:)]; %equacao especifica para esse trafo

figure
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on
plot(xx*1e6,yy(1,:),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx*1e6,yy(2,:),'Color',jetcustom(2,:),'LineWidth',1.5)
plot(xx*1e6,yy(3,:),'Color',jetcustom(3,:),'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
legend({'$i_{Ldab_a}$','$i_{Ldab_b}$','$i_{Ldab_c}$'},'Location','southeast','FontSize', 16)






%% tensao e corrente sobre o indutor


fVLdab = matlabFunction(dwb_L, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

yy = fVLdab(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy,yy(1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

yy2 = fidab(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = [yy2,yy2(1)];


fid = fopen('C1Trace00002.dat');
format long
cell_data= textscan(fid,'%f%f','Delimiter',' ','headerLines',1);
my_data1 = cat(2,cell_data{:});
fclose(fid);

fid = fopen('C2Trace00000.dat');
format long
cell_data= textscan(fid,'%f%f','Delimiter',' ','headerLines',1);
my_data2 = cat(2,cell_data{:});
fclose(fid);

figure

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
L11 = plot(nan, nan,'-','LineWidth',1.5,'Color',[0 0 0]);
L22 = plot(nan, nan,':','LineWidth',1.5,'Color',[0 0 0]);

stairs(xx*1e6,yy,':','Color',jetcustom(1,:), 'LineWidth',1.5)
plot((-my_data2(1,1)+my_data2(:,1))*1e6,my_data2(:,2),'-','Color',jetcustom(1,:), 'LineWidth',1.5)
ylabel('$V_{La}\,$[V]')
yyaxis right
plot(xx*1e6,yy2,':','Color',jetcustom(2,:), 'LineWidth',1.5)
plot((-my_data1(1,1)+my_data1(:,1))*1e6,my_data1(:,2),'-','Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
% xlim([0 1/fs_num])
xlabel('$t[\mu$s]')
ylabel('$i_{La}\,$[A]')
legend({'Experimental','Theoric'},'Location','best','FontSize', 14)
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
f_save_figure(append('figure\comp\fig1.pdf'))


%% corrente primario e secundario

fid = fopen('C1Trace00005.dat');
format long
cell_data= textscan(fid,'%f%f','Delimiter',' ','headerLines',1);
my_data1 = cat(2,cell_data{:});
fclose(fid);

fid = fopen('C3Trace00004.dat');
format long
cell_data= textscan(fid,'%f%f','Delimiter',' ','headerLines',1);
my_data2 = cat(2,cell_data{:});
fclose(fid);

fil= matlabFunction(x0s(1,:), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiL = matlabFunction(x0s(4,:), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

yy = fil(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = fiL(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);


figure


cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
L11 = plot(nan, nan,'-','LineWidth',1.5,'Color',[0 0 0]);
L22 = plot(nan, nan,':','LineWidth',1.5,'Color',[0 0 0]);

plot(xx*1e6,yy,':','Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx*1e6,-yy2,':','Color',jetcustom(2,:),'LineWidth',1.5)

kk = 10;
text(xx(kk)*1e6,yy(kk)+2,'$i_p$','Interpreter', 'Latex','FontSize', 16) 
kk = 5;
text(xx(kk)*1e6,-yy2(kk)+2,'$-i_s$','Interpreter', 'Latex','FontSize', 16) 

plot((-my_data1(1,1)+my_data1(:,1))*1e6,my_data1(:,2),'-','Color',jetcustom(1,:), 'LineWidth',1.5)
plot((-my_data2(1,1)+my_data2(:,1))*1e6,my_data2(:,2),'-','Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
legend({'Experimental','Theoric'},'Location','best','FontSize', 14)
ylim([-15 15])
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
f_save_figure(append('figure\comp\fig2.pdf'))

%% correntes nos indutores

fid = fopen('C1Trace00006.dat');
format long
cell_data= textscan(fid,'%f%f','Delimiter',' ','headerLines',1);
my_data1 = cat(2,cell_data{:});
fclose(fid);

fid = fopen('C3Trace00005.dat');
format long
cell_data= textscan(fid,'%f%f','Delimiter',' ','headerLines',1);
my_data2 = cat(2,cell_data{:});
fclose(fid);

fid = fopen('C4Trace00002.dat');
format long
cell_data= textscan(fid,'%f%f','Delimiter',' ','headerLines',1);
my_data3 = cat(2,cell_data{:});
fclose(fid);

fl= matlabFunction(x0s(1:3,:), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
yy = fl(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);


figure
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on

L11 = plot(nan, nan,'-','LineWidth',2.5,'Color',[0 0 0]);
L22 = plot(nan, nan,'--.','LineWidth',1.5,'Color',[0 0 0]);

plot(xx*1e6,yy(1,:),':','Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx*1e6,yy(2,:),':','Color',jetcustom(2,:),'LineWidth',1.5)
plot(xx*1e6,yy(3,:),':','Color',jetcustom(3,:),'LineWidth',1.5)


plot((-my_data1(1,1)+my_data1(:,1))*1e6,my_data1(:,2),'-','Color',jetcustom(1,:), 'LineWidth',1.5)
plot((-my_data2(1,1)+my_data2(:,1))*1e6,my_data2(:,2),'-','Color',jetcustom(3,:), 'LineWidth',1.5)
plot((-my_data3(1,1)+my_data3(:,1))*1e6,my_data3(:,2),'-','Color',jetcustom(2,:), 'LineWidth',1.5)

hold off
grid on
grid minor
set(gca, 'FontSize', 20)
legend({'Experimental','Theoric'},'Location','best','FontSize', 14)
f_save_figure(append('figure\comp\fig3.pdf'))
