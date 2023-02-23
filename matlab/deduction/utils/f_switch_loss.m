function [ZVS_p, ZVS_s, prim_loss, sec_loss] = f_switch_loss(I_sw_p_on,I_sw_s_on,Vp,d,fs)

    load('f_fitted_off.mat')
    load('f_fitted_on.mat')

    if (I_sw_p_on < 0)
        ZVS_p = 1;
    else
        ZVS_p = 0;
    end

    if (I_sw_s_on < 0)
        ZVS_s = 1;
    else
        ZVS_s = 0;
    end

    if (ZVS_p == 1)
        e_p_turn_on = 0;
        e_p_turn_off = f_fitted_off(-I_sw_p_on,Vp);
    else
        e_p_turn_on = f_fitted_on(I_sw_p_on,Vp);
        e_p_turn_off = 0;
    end
    
    if (ZVS_s == 1)
        e_s_turn_on = 0;
        e_s_turn_off = f_fitted_off(-I_sw_s_on,Vp*d);
    else
        e_s_turn_on = f_fitted_on(I_sw_s_on,Vp*d);
        e_s_turn_off = 0;
    end

    prim_loss = fs*(e_p_turn_on+e_p_turn_off)*6;
    sec_loss = fs*(e_s_turn_on+e_s_turn_off)*6;
end