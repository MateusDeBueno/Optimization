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
phi_num = deg2rad(0.01);

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
    x0s = simplify(Tclxx\subs(xcl, x0, x0x));
    dx0s = Tclxx\B*scl; %derivadas dos estados, indutor e trafo secundario
        
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
    
    %% algumas coisas pra plotar
    for ii=1:length(x0s)-1
        vSS(:,ii) = subs(vS,dx,dx0s(:,ii)); %tensao aplicada na bobinas do sec da fase A
        vPP(:,ii) = subs(vP,dx,dx0s(:,ii)); %tensao aplicada na bobinas do prim da fase A
        vLLdab(:,ii) = subs(vLdab,dx,dx0s(:,ii)); %tensao aplicada no indutor da fase A
        HBB(:,ii) = subs(iHB,x,x0s(:,ii)); %correntes nos 3 half bridge primario
        hbb(:,ii) = subs(ihb,x,x0s(:,ii)); %correntes nos 3 half bridge secundario
        idab_plot(:,ii) = subs(ild,x,x0s(:,ii));
    end 

    vSS = [vSS,vSS(:,1)];
    vPP = [vPP,vPP(:,1)];
    vLLdab = [vLLdab,vLLdab(:,1)];
    HBB = [HBB,HBB(:,1)];
    hbb = [hbb,hbb(:,1)];
    idab_plot = [idab_plot,idab_plot(:,1)];

    vSS = simplify(vSS);
    vPP = simplify(vPP);
    vLLdab = simplify(vLLdab);
    HBB = simplify(HBB);
    hbb = simplify(hbb);
    idab_plot = simplify(idab_plot);
    
    sec_switch_plot = ones(1,length(ts))*sec_switch;
    sf_plot = [sf,sf(:,1)];

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
    
    dB_tr = vSS(1,:)/(N_tr*Ac_tr);
    
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


    
    f_plot = matlabFunction([vPP;vSS;vLLdab;hbb;HBB;idab_plot;x0s;ts;sec_switch_plot;sf_plot],'vars',{L1,L2,Ldab,M,Vi,d,fs,phi});
    f = matlabFunction([hbrm;HBrm;Ip;Is;iiRMS;iiME;ioRMS;ioME;Pm;idrm;ilrm;iLrm;iSwPrm;iSwSrm;P_core_tr;Bpk_tr;P_core_L;Bpk_L],...
        'vars',{L1,L2,Ldab,M,Vi,d,fs,phi,ki_tr,b_tr,a_tr,Ac_tr,N_tr,Ve_tr,ki_L,b_L,a_L,Ac_L,N_L,Ve_L});
%     fh = matlabFunction([ilrm_cn;iLrm_cn;idrm_cn],'vars',{L1,L2,Ldab,M,Vi,d,fs,phi,nn});

% ssssdgbgv

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
fidab_plot = matlabFunction(idab_plot, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%% indutor
yy = fidab_plot(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = fvLLdab(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

fase = 3;

figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
plot(xx*1e6,yy(fase,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
ylabel('$V_{la}\,$[V]')
yyaxis right
stairs(xx*1e6,yy2(fase,:),'Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i_{la}\,$[A]')
yyaxis left
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
title('Indutor')


%% primario
yy = fx0s(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = fvPP(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

fase = 3;

figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
plot(xx*1e6,yy(fase,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
ylabel('$V_{la}\,$[V]')
yyaxis right
stairs(xx*1e6,yy2(fase,:),'Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i_{la}\,$[A]')
yyaxis left
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
title('primario')



%% secundario
yy = fx0s(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
yy2 = fvPP(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);
xx = fts(L1_num,L2_num,Ldab_num,M_num,Vi_num,d_num,fs_num,phi_num);

fase = 3;

figure
cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on
plot(xx*1e6,yy(3+fase,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
ylabel('$V_{la}\,$[V]')
yyaxis right
stairs(xx*1e6,yy2(fase,:),'Color',jetcustom(2,:), 'LineWidth',1.5)
hold off
grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i_{la}\,$[A]')
yyaxis left
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
set(gca, 'FontSize', 20)
xlabel('$t[\mu$s]')
ylabel('$i\,$[A]')
title('secundario')