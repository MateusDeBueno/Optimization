function [output,output_harmonic,output_phi] = dabYD(phi_num)
    syms phi fs Vi Vo Ldab Ld1 Ld2 Lm n dt Po k real    
    
    Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
    Tclm = Tcl(1:2,:); %ignorar nivel zero
    Tclx = kron(eye(2), Tclm);
    
    Ts = 1/fs;
    
    [sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi);
    
    Tds = [1 -1 0;0 1 -1;-1 0 1];
    sf_line = Tds*sf_s;
    
    sf = [sf_p; sf_line]; %a1b1c1 to a2b2c2

    scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
    ts = simplify(ang)*Ts/(2*pi);
    tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos
        
    %% Resolvendo circuito monofasico equivalente [MUDAR]
    L = [Ldab+Ld1; Ld2/n^2];
    Udc = [Vi; Vo];
    u = Udc.*[1; 1/n];
    
    x = sym('x_', [2,1], 'real');
    dx = sym('dx_', [2,1], 'real'); %dx1 -> iLdab ou iLl1 (em serie) dx2 -> iLl2
    eq = sym('eq', [2,1], 'real');
    s = sym('s', [2,1], 'real');
    
    eq(1) = u(1)*s(1) == L(1)*dx(1) + Lm*(dx(1)-dx(2)*n);
    eq(2) = u(2)*s(2) == Lm*(dx(1)-dx(2)*n) - L(2)*dx(2)*n;
    
    dxs = simplify(struct2array(solve(eq, dx))).';
    M = equationsToMatrix(dxs, [x; s]);
    A = M(:,1:2);
    B = M(:,3:4); %derivadas de iLdab e iLd
    
    A = kron(A, eye(2));
    B = kron(B, eye(2));
    
    %% Obter valores de regime permanente
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
    
    %% Corrente eficaz nos estados
    [x_rms,~] = rms_and_mean(derivadas,x0s,ts,(1:12),(1:12));
    IL_rms = x_rms(1);
    Itrfsec_rms = x_rms(4);

    %% Fourier dos estados
    IL_rms_c_k = ck_fourier(ts,derivadas(1,:),x0s(1,:)); %DESCOMENTAR
    Itrf_sec_rms_c_k = ck_fourier(ts,derivadas(4,:),x0s(4,:)); %DESCOMENTAR
    
    %% Corrente do half bridge do primario [MUDAR]
    derivadas_hb_p = derivadas(1:3,:);
    pts_inics_hb_p = x0s(1:3,:);
    
    %% Corrente do half bridge do secundario [MUDAR]
    derivadas_hb_s = derivadas([4 5 6],:) - derivadas([6 4 5],:);
    pts_inics_hb_s = x0s([4 5 6],:) - x0s([6 4 5],:);

    %% Corrente de entrada
    target = [1;0;0];
    [etapas] = pega_etapa(sf_p,target);
    [Iin_rms,Iin_med] = rms_and_mean(derivadas_hb_p(1,:),pts_inics_hb_p(1,:),ts,etapas,etapas);
    
    %% Corrente de saida
    target = [1;0;0];
    [etapas] = pega_etapa(sf_s,target);
    [Iout_rms,Iout_med] = rms_and_mean(derivadas_hb_s(1,:),pts_inics_hb_s(1,:),ts,etapas,etapas);
    
    power_equation = Iout_med==Po/Vo;
    [phi_sol,~,phi_cond] = solve(power_equation,phi,'ReturnConditions',true);

    if (length(phi_sol)==1)
        output_phi = [phi_sol;phi_sol;phi_cond;phi_cond];
    else
        output_phi = [phi_sol;phi_cond];
    end
    %% Corrente nas chaves
    [Isw_p_rms,~] = rms_and_mean(derivadas_hb_p(1,:),pts_inics_hb_p(1,:),ts,(1:6),(1:12));
    
    [Isw_s_rms,~] = rms_and_mean(derivadas_hb_s(1,:),pts_inics_hb_s(1,:),ts,(1:6),(1:12));
    
    %% corrente de comutacao
    %primario
    Ip = pts_inics_hb_p(1,1);
    
    %secundario
    Is = -pts_inics_hb_s(1,sec_switch);

    %% output
    output = [IL_rms,Itrfsec_rms,Iin_med,Iin_rms,Iout_med,Iout_rms,Isw_p_rms,Isw_s_rms,Ip,Is];
    output_harmonic = [IL_rms_c_k,Itrf_sec_rms_c_k];
end