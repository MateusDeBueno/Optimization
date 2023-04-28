function [f,fh,f_plot] = full_DfD(phi_num)
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
    iHB = [iL(3); iL(1:2)]-iL;
    diHB = [diL(3); diL(1:2)]-diL;
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
    ilrm_cn = ck_fourier(ts,dx0s(1,:),x0s(1,:)); %DESCOMENTAR
    iLrm_cn = ck_fourier(ts,dx0s(4,:),x0s(4,:)); %DESCOMENTAR  
    
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
    idrm_cn = ck_fourier(ts,didab(1,:),idab(1,:)); %DESCOMENTAR
    
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
    fh = matlabFunction([ilrm_cn;iLrm_cn;idrm_cn],'vars',{L1,L2,Ldab,M,Vi,d,fs,phi,nn});
end