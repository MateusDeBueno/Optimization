function [cond_loss] = f_cond_loss(Isw_rms_p,Isw_rms_s)   
    
    Ron_p = 80e-3;
    Ron_s = Ron_p;

    p_primary = Ron_p*Isw_rms_p^2;
    p_secundary = Ron_s*Isw_rms_s^2;

    cond_loss = (p_primary+p_secundary)*6;
end

