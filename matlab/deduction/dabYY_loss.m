function [Pt] = dabYY_loss(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo,dab)

    if (phi>0 && phi<pi/3)
        Isw_p = dab.YYmenor60.f_Ip(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
        Isw_s = dab.YYmenor60.f_Is(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
        Isw_p_rms = dab.YYmenor60.f_Isw_p_rms(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
        Isw_s_rms = dab.YYmenor60.f_Isw_s_rms(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
    elseif (phi < pi/2)
        Isw_p = dab.YYmaior60.f_Ip(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
        Isw_s = dab.YYmaior60.f_Is(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
        Isw_p_rms = dab.YYmaior60.f_Isw_p_rms(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
        Isw_s_rms = dab.YYmaior60.f_Isw_s_rms(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);
    end

    %% Perdas de comutacao
    Psw_p = sw_loss(fs, Vi, Isw_p, dab);
    Psw_s = sw_loss(fs, Vo, Isw_s, dab);
    
    %% Perdas de conducao
    Pcond_p = dab.sw.Rds_on*Isw_p_rms^2;
    Pcond_s = dab.sw.Rds_on*Isw_s_rms^2;
    
    %% Perdas no cobre do indutor
    Pcu_trf = 0;

    %% Perdas no nucleo do indutor
    Pco_trf = 0;

    %% Perdas no cobre do transformador
    Pcu_ind = 0;

    %% Perdas no nucleo do transformador
    Pco_ind = 0;

    %% Perdas totais
    Pt = (Psw_p + Psw_s + Pcond_p + Pcond_s)*6 + Pcu_trf + Pco_trf + (Pcu_ind + Pco_ind)*3;
end