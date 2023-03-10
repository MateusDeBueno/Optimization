function [x0s, ts, idab, hb, HB, Ip, Is, iME, idrm, ilrm, iLrm, iSwPrm, iSwSrm,B] = simplify_DfY(phi_num)
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
    
    %corrente no Ldab
    ild = il - [il(3); il(1:2)];
    dild = dil - [dil(3); dil(1:2)];
    
    %definicao das tensoes do indutor acoplado
    vP = + L1*dil + M*diL;
    vS = + M*dil + L2*diL;
    
    %definicao das tensoes nos indutores do dab
    vLdab = Ldab*dild;
    
    Td = [1 -1 0; 0 1 -1; -1 0 1];
    
    %malhas
    m_p = Td*u(1:3) == Td*vLdab + vP; %malha primario
    m_s = Td*u(4:6) == Td*vS; %malha secundario
    
    %usa 4 malhas
    eq(1:2) = m_p(1:2);
    eq(3:4) = m_s(1:2);
    
    %% completar ilc e iLC
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
    
%     Axx = simplify(Tclxz*As*inv(Tclxz));
%     Bxx = simplify(Tclxz*Bs*inv(Tclxz)); %a pinv deixa uns numeros feios
%     Ax = [Axx(:,1:2),Axx(:,4:5)];
%     Bx = [Bxx(:,1:2),Bxx(:,4:5)];
%     A = [Ax(1:2,:);Ax(4:5,:)];
%     B = [Bx(1:2,:);Bx(4:5,:)];
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
    
    %% 
    %obtido os valores de corrente na bobina do primario e secundario
    %eh obtido todas as outras corrente
    
    %corrente no Ldab
    for ii=1:length(x0s)-1
        idab(:,ii) = subs(ild(1),x(1:3),x0s(1:3,ii));
        didab(:,ii) = subs(dild(1),dx(1:3),dx0s(1:3,ii));
    end
    
    hb = idab;
    d_hb = didab;
    HB = -x0s(4:6,:);
    d_HB = -dx0s(4:6,:);

    %% Corrente Ldab
    [idrm,~] = rms_and_mean(d_hb(1,:),hb(1,:),ts,[1:12],[1:12]);
    
    %% Corrente nos estados
    [ilrm,~] = rms_and_mean(dx0s(1,:),x0s(1,:),ts,[1:12],[1:12]);
    [iLrm,~] = rms_and_mean(dx0s(4,:),x0s(4,:),ts,1:12,1:12);

    %% Corrente de entrada
%     target = [1;0;0];
%     [etapas] = pega_etapa(sf_p,target);
%     [irm,ime] = rms_and_mean(d_hb(1,:),hb(1,:),ts,etapas,etapas);
    
    %% Corrente de saida
    target = [1;0;0];
    [etapas] = pega_etapa(sf_s,target);
    [~,iME] = rms_and_mean(d_HB(1,:),HB(1,:),ts,etapas,etapas);
    
    %% Corrente nas chaves
    [iSwPrm,~] = rms_and_mean(d_hb(1,:),hb(1,:),ts,1:6,1:12);
    [iSwSrm,~] = rms_and_mean(d_HB(1,:),HB(1,:),ts,1:6,1:12);  

    %% corrente de comutacao
    Ip = hb(1,1);
    Is = -HB(1,sec_switch);

end