function Psw = sw_loss(fs, Vblock, Isw, dab)
    
    if Isw < 0
        Esw_on = 0;
    else
        Esw_on = dab.sw.fiton.f_fitted_on(Isw,Vblock);
    end
    
    %  Energia da Saída de Condução
    if Isw > 0
        Esw_off = 0;
    else
        Esw_off = dab.sw.fitoff.f_fitted_off(-Isw,Vblock);
    end
    
    Psw = fs*(Esw_on + Esw_off);
end