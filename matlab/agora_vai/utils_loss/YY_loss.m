function [efc,ptot,Pm,cSw_p,cSw_s,sSw_p,sSw_s,Pil_cu,PiL_cu,PiLd_cu,Psc,P_core_L,P_core_tr] = YY_loss(l)
    %% equations
    out = f_equationsYY(l);
    C = num2cell(out);
    [hbrm,HBrm,Ip,Is,Pm,~,~,~,iSwPrm,iSwSrm,P_core_tr,~,P_core_L,~] = C{:};
    

%     if (Ip>0)
%         l.pr.phi = l.pr.phi - (l.pr.dt/2)*l.pr.fs*2*pi;
%         out = f_equationsYY(l);
%         C = num2cell(out);
%         [hbrm,HBrm,Ip,Is,Pm,~,~,~,iSwPrm,iSwSrm,P_core_tr,~,P_core_L,~] = C{:};
%     end
%     if (Is>0)
%         l.pr.phi = l.pr.phi + (l.pr.dt/2)*l.pr.fs/(2*pi);
%         out = f_equationsYY(l);
%         C = num2cell(out);
%         [hbrm,HBrm,Ip,Is,Pm,~,~,~,iSwPrm,iSwSrm,P_core_tr,~,P_core_L, ] = C{:};
%     end


    %% semi conductor loss
    
    cSw_p = 6*l.sw.Ronp*iSwPrm^2;
    cSw_s = 6*l.sw.Rons*iSwSrm^2;
    
    sSw_p = 6*f_sw_loss(l.pr.fs, l.pr.Vi, Ip, l);
    sSw_s = 6*f_sw_loss(l.pr.fs, l.pr.Vi*l.pr.d, Is, l);
    
    %% wire loss
    Pil_cu = 0;
    PiL_cu = 0;
    PiLd_cu = 0;
    for nn=1:100
        [R_ac_p, R_ac_s] = f_get_resistance_trf(nn,0.8,l);
        RacL = f_get_resistance_ind(nn,0.8,l);
        
        out = fh_equationsYY(l,nn);
        C = num2cell(out);
        [ilrm_cn,iLrm_cn,idrm_cn] = C{:};
        
        Pil_cu = Pil_cu + R_ac_p*ilrm_cn^2;
        PiL_cu = PiL_cu + R_ac_s*iLrm_cn^2;
        PiLd_cu = PiLd_cu + RacL*idrm_cn^2;
    end
    
    Pil_cu = 3*Pil_cu;
    PiL_cu = 3*PiL_cu;
    PiLd_cu = 3*PiLd_cu;
    
    %% serie capacitor
    Psc = 3*l.sC.Rac100*hbrm^2 + 3*l.sC.Rac100*HBrm^2;
    
    %%
    
    ptot = cSw_p+cSw_s+sSw_p+sSw_s+Pil_cu+PiL_cu+PiLd_cu+P_core_L+P_core_tr+Psc;
    
    efc = (Pm-ptot)/Pm;
end