function [out] = fDinY_harm(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num)
    
    if (-pi <= phi_num && phi_num <= -2*pi/3)
        out = trafo.DinY.str_h.f1(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num);
    elseif (-2*pi/3 < phi_num && phi_num <= -pi/3)
        out = trafo.DinY.str_h.f2(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num);
    elseif (-pi/3 < phi_num && phi_num <= 0)
        out = trafo.DinY.str_h.f3(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num);
    elseif (0 < phi_num && phi_num <= pi/3)
        out = trafo.DinY.str_h.f4(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num);
    elseif (pi/3 < phi_num && phi_num <= 2*pi/3)
        out = trafo.DinY.str_h.f5(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num);
    elseif (2*pi/3 < phi_num && phi_num <= pi)
        out = trafo.DinY.str_h.f6(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, phi_num, fs_num, Vi_num, Vo_num, k_num);
    end
end