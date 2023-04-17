clear; close all; clc;
%%

trafo = 'oDY';

color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue
color_blue_eth = [33, 92, 175]/255;

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

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61.5e-6;
Ld1_num = 2e-6;
n_num = 1;
Ld2_num = 2e-6;
Lm_num = .5e-3;
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);
phi_num = deg2rad(-15);

%%

%cria intervalo de angulo
% intervalo = [0 deg2rad(60)];
intervalo = [deg2rad(-60) 0];

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
    
    Tclxx = kron(eye(2), Tcl);

    A = Tclxx*As/Tclxx;
    B = Tclxx*Bs/Tclxx;
    
    %% obter funcao de comutacao, depende de phi
    [sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi);
    sf = sf+0.5;
    scl = kron(eye(2), Tcl)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
    ts = simplify(ang)*Ts/(2*pi);
    tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos
    
    %% Obter valores de regime permanente
    g = simplify(expm(A*dt));
    h = B*dt;
    xcl = sym('x_', [6,length(ts)], 'real');
    x0 = sym('x0_',[6,1], 'real');
    
    xcl(:,1) = x0;
    for i=1:length(ts)-1
        dts = (ts(i+1)-ts(i));
        xcl(:,i+1) = simplify(subs(g, dt, dts)*xcl(:,i) + subs(h, dt, dts)*scl(:,i));
    end
    
    x0x = struct2array(solve(xcl(:,1) == -xcl(:,7), x0)).';
    x0s = simplify(inv(Tclxx)*subs(xcl, x0, x0x));
    dx0s = inv(Tclxx)*B*scl; %derivadas dos estados, indutor e trafo secundario
        
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

% fP_core_tr(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.tr.N,l.tr.Ac,l.tr.Ve,l.tr.a,l.tr.b,l.tr.ki)
% fBpk_tr(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.tr.N,l.tr.Ac)
% fP_core_L(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.L.N,l.L.Ac,l.L.Ve,l.L.a,l.L.b,l.L.ki)
% fBpk_L(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.L.N,l.L.Ac)

% Pse = l.tr.Ve*l.tr.kc*(l.pr.fs)^(l.tr.a)*fBpk_tr(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num,l.tr.N,l.tr.Ac)^(l.tr.b);


%%

for ii=1:length(x0s)-1
    vSS(ii) = subs(vS(1),dx,dx0s(:,ii));
    vPP(ii) = subs(vP(1),dx,dx0s(:,ii));
    vLLdab(ii) = subs(vLdab(1),dx,dx0s(:,ii));
    HBB(:,ii) = subs(iHB,x,x0s(:,ii));
    hbb(:,ii) = subs(ihb,x,x0s(:,ii));
end


%% criar uma porrada de funcoes
syms Po real

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

fvSS = matlabFunction(vSS, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fvPP = matlabFunction(vPP, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fvLLdab = matlabFunction(vLLdab, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fHBB = matlabFunction(HBB, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fhbb = matlabFunction(hbb, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%% tensao e corrente sobre o indutor

yy = fvLLdab(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy,yy(1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = fidab(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = [yy2,yy2(1)];

figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
stairs(xx*1e6,yy,'Color',jetcustom(1,:), 'LineWidth',1.5)
ylabel('$V_{la}\,$[V]')
yyaxis right
plot(xx*1e6,yy2,'Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i_{la}\,$[A]')
ylim([-5 5])
yyaxis left
ylim([-300 300])
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
file_name = append('figure\finalCap2\VLdabiLdab_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% tensao e corrente no primario



yy = fvPP(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy,yy(1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = fx0s(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = yy2(1,:);
figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
stairs(xx*1e6,yy,'Color',jetcustom(1,:), 'LineWidth',1.5)
ylabel('$V_{L1_{a}}\,$[V]')
yyaxis right
plot(xx*1e6,yy2,'Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i_{tra}\,$[A]')
ylim([-3 3])
yyaxis left
ylim([-300 300])
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
file_name = append('figure\finalCap2\Va_ia_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% tensao e corrente no secundario

yy = fvSS(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy,yy(1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = fx0s(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = yy2(4,:);
figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
stairs(xx*1e6,yy,'Color',jetcustom(1,:), 'LineWidth',1.5)
ylabel('$V_{L2_{A}}\,$[V]')
yyaxis right
plot(xx*1e6,yy2,'Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i_{trA}\,$[A]')
ylim([-3 3])
yyaxis left
ylim([-300 300])
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
file_name = append('figure\finalCap2\VAiA_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% corrente no primario e secundario

yy = fx0s(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on
plot(xx*1e6,yy(1,:),'Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx*1e6,-yy(4,:), 'Color',jetcustom(2,:),'LineWidth',1.5)
hold off
ylim([-3 3])
grid on
grid minor
legend({'$i_{la}$','$-i_{La}$'},'Location','best','FontSize', 16)
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
file_name = append('figure\finalCap2\iliL_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% corrente de saida

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

target = [1;0;0];
[etapas] = pega_etapa(sf_s,target);
[ioRMS,ioME] = rms_and_mean(dHB(1,:),HB(1,:),ts,etapas,etapas);

target_1n = [0;1;1]; %1o p neg
target_2n = [1;0;1]; %2o p neg
target_3n = [1;1;0]; %3o p neg
target_1p = [1;0;0]; %1o p pos
target_2p = [0;1;0]; %2o p pos
target_3p = [0;0;1]; %3o p pos

[etapas_1p] = pega_etapa(sf_s,target_1p);
[etapas_2p] = pega_etapa(sf_s,target_2p);
[etapas_3p] = pega_etapa(sf_s,target_3p);
[etapas_1n] = pega_etapa(sf_s,target_1n);
[etapas_2n] = pega_etapa(sf_s,target_2n);
[etapas_3n] = pega_etapa(sf_s,target_3n);

%ajeita etapa
etapas_1p = [etapas_1p max(etapas_1p)+1];
etapas_2p = [etapas_2p max(etapas_2p)+1];
etapas_3p = [etapas_3p max(etapas_3p)+1];
etapas_1n = [etapas_1n max(etapas_1n)+1];
etapas_2n = [etapas_2n max(etapas_2n)+1];
etapas_3n = [etapas_3n max(etapas_3n)+1];



yy = fHBB(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy, yy(:,1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
xx = [xx, xx(end)+xx(2:end)];

i_adapt = yy(1,etapas_1p);
i_vert = [i_adapt(end), i_adapt(1)];
%%
linewid = 2.2;

figure
hold on

plot(nan,nan,'-.','Color',jetcustom(1,:),'LineWidth',1.5)
plot(nan,nan,'-.','Color',jetcustom(2,:),'LineWidth',1.5)
plot(nan,nan,'-.','Color',jetcustom(3,:),'LineWidth',1.5)

for kk=2:2:12
    plot([xx(kk), xx(kk)]*1e6,i_vert,'Color',color_blue_eth,'LineWidth',linewid)
end

plot(xx(1,etapas_1p)*1e6,yy(1,etapas_1p),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_2p)*1e6,yy(2,etapas_2p),'Color',color_blue_eth,'LineWidth',linewid)
etapas_3p1 = etapas_3p(2:3);
etapas_3p2 = [etapas_3p(1) etapas_3p(1)+1];
plot(xx(1,etapas_3p1)*1e6,yy(3,etapas_3p1),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_3p2)*1e6,yy(3,etapas_3p2),'Color',color_blue_eth,'LineWidth',linewid)

plot(xx(1,etapas_1n)*1e6,-yy(1,etapas_1n),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_2n)*1e6,-yy(2,etapas_2n),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_3n)*1e6,-yy(3,etapas_3n),'Color',color_blue_eth,'LineWidth',linewid)


plot(xx(1:13)*1e6,-yy(1,:),'-.','Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx(1:13)*1e6,-yy(2,:),'-.','Color',jetcustom(2,:),'LineWidth',1.5)
plot(xx(1:13)*1e6,-yy(3,:),'-.','Color',jetcustom(3,:),'LineWidth',1.5)

x_points = [xx(2), xx(2), xx(4), xx(4)]*1e6;  
y_points = [-10, 10, 10, -10];
color = [0, 0, 1];

a = fill(x_points, y_points, color, 'LineStyle', 'none');
a.FaceAlpha = 0.1;

text(xx(4)*1e6,3.5,'$\leftarrow 2^{\circ}$ and $3^{\circ}$ stages','Interpreter', 'Latex','FontSize', 16) 


hold off
ylim([-4 4])
grid on
grid minor
legend({'$-i_{A}$','$-i_{B}$','$-i_{C}$','$i_{out}$'},'Location','southeast','FontSize', 16)
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')

file_name = append('figure\finalCap2\iout_part1_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');
%%
figure
hold on

plot(nan,nan,'-.','Color',jetcustom(1,:),'LineWidth',1.5)
plot(nan,nan,'-.','Color',jetcustom(2,:),'LineWidth',1.5)
plot(nan,nan,'-.','Color',jetcustom(3,:),'LineWidth',1.5)

for kk=2:2:12
    plot([xx(kk), xx(kk)]*1e6,i_vert,'Color',color_blue_eth,'LineWidth',linewid)
end

plot(xx(1,etapas_1p)*1e6,yy(1,etapas_1p),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_2p)*1e6,yy(2,etapas_2p),'Color',color_blue_eth,'LineWidth',linewid)
etapas_3p1 = etapas_3p(2:3);
etapas_3p2 = [etapas_3p(1) etapas_3p(1)+1];
plot(xx(1,etapas_3p1)*1e6,yy(3,etapas_3p1),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_3p2)*1e6,yy(3,etapas_3p2),'Color',color_blue_eth,'LineWidth',linewid)

plot(xx(1,etapas_1n)*1e6,-yy(1,etapas_1n),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_2n)*1e6,-yy(2,etapas_2n),'Color',color_blue_eth,'LineWidth',linewid)
plot(xx(1,etapas_3n)*1e6,-yy(3,etapas_3n),'Color',color_blue_eth,'LineWidth',linewid)


plot(xx(1:13)*1e6,yy(1,:),'-.','Color',jetcustom(1,:),'LineWidth',1.5)
plot(xx(1:13)*1e6,yy(2,:),'-.','Color',jetcustom(2,:),'LineWidth',1.5)
plot(xx(1:13)*1e6,yy(3,:),'-.','Color',jetcustom(3,:),'LineWidth',1.5)


hold off
ylim([-4 4])
grid on
grid minor
legend({'$i_{A}$','$i_{B}$','$i_{C}$','$i_{out}$'},'Location','southeast','FontSize', 16)
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')

file_name = append('figure\finalCap2\iout_part2_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% corrente na chave primario
[iSwPrm,~] = rms_and_mean(dhb(1,:),hb(1,:),ts,1:6,1:12);
[iSwSrm,~] = rms_and_mean(dHB(1,:),HB(1,:),ts,1:6,1:12);  


yy = fhbb(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy, yy(:,1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

xxp1 = [xx(1:7), xx(7)];
yyp1 = [yy(1,1:7), 0];
xxp2 = [xx(7) xx(end)];
yyp2 = [0 0];

xxp = [xxp1 xxp2];
yyp = [yyp1 yyp2];

figure
hold on
plot(xxp*1e6,yyp,'Color',color2,'LineWidth',2)
plot(xx*1e6,yy(1,:),'-.','Color',color1,'LineWidth',1.5)
hold off
ylim([-5 5])
grid on
grid minor
legend({'$i_{sw-p}$','$i_{a}$'},'Location','best','FontSize', 16)
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
file_name = append('figure\finalCap2\isw_p_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% corrente na chave secundario


yy = fHBB(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy, yy(:,1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

xxp0 = [xx(1) xx(2)];
yyp0 = [0 0];
xxp1 = [xx(sec_switch:sec_switch+6), xx(sec_switch+6)];
yyp1 = [yy(1,sec_switch:sec_switch+6), 0];
xxp2 = [xx(sec_switch+6) xx(end)];
yyp2 = [0 0];


xxp = [xxp0 xxp1 xxp2];
yyp = [yyp0 yyp1 yyp2];

figure
hold on
plot(xxp*1e6,-yyp,'Color',color2,'LineWidth',2)
plot(xx*1e6,-yy(1,:),'-.','Color',color1,'LineWidth',1.5)
hold off
ylim([-3 3])
grid on
grid minor
legend({'$i_{sw-s}$','$-i_{A}$'},'Location','best','FontSize', 16)
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
file_name = append('figure\finalCap2\isw_s_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%%
latex(ioME)
latex(Ip)
latex(Is)



%%
figure
sf_p


[iSwPrm,~] = rms_and_mean(dhb(1,:),hb(1,:),ts,1:6,1:12);
[iSwSrm,~] = rms_and_mean(dHB(1,:),HB(1,:),ts,1:6,1:12);  


yy = fhbb(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy = [yy, yy(:,1)];
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

teste = yy.*[sf_p, sf_p(:,1)];


figure
plot(xx, teste)

% 
% xxp1 = [xx(1:7), xx(7)];
% yyp1 = [yy(1,1:7), 0];
% xxp2 = [xx(7) xx(end)];
% yyp2 = [0 0];
% 
% xxp = [xxp1 xxp2];
% yyp = [yyp1 yyp2];
% 
% figure
% hold on
% plot(xxp*1e6,yyp,'Color',color2,'LineWidth',2)
% plot(xx*1e6,yy(1,:),'-.','Color',color1,'LineWidth',1.5)
% hold off
% ylim([-5 5])
% grid on
% grid minor
% legend({'$i_{sw-p}$','$i_{a}$'},'Location','best','FontSize', 16)
% set(gca, 'FontSize', 20)
% xlabel('$t[\mu$s]')
% ylabel('$i\,$[A]')
% file_name = append('figure\finalCap2\isw_p_',trafo,'.pdf');
% exportgraphics(gca,file_name,'ContentType','vector');
