function out_phi = fDinY_phi(trafo,Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num)
    
    if (-pi <= phi_num && phi_num <= -2*pi/3)
        eqx = trafo.DinY.fphi.f1(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num);
        if(eqx(3))
            out_phi = eqx(1);
        elseif(eqx(4))
            out_phi = eqx(2);
        end
    elseif (-2*pi/3 < phi_num && phi_num <= -pi/3)
        eqx = trafo.DinY.fphi.f2(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num);
        if(eqx(3))
            out_phi = eqx(1);
        elseif(eqx(4))
            out_phi = eqx(2);
        end
    elseif (-pi/3 < phi_num && phi_num <= 0)
        eqx = trafo.DinY.fphi.f3(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num);
        if(eqx(3))
            out_phi = eqx(1);
        elseif(eqx(4))
            out_phi = eqx(2);
        end
    elseif (0 < phi_num && phi_num <= pi/3)
        eqx = trafo.DinY.fphi.f4(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num);
        if(eqx(3))
            out_phi = eqx(1);
        elseif(eqx(4))
            out_phi = eqx(2);
        end
    elseif (pi/3 < phi_num && phi_num <= 2*pi/3)
        eqx = trafo.DinY.fphi.f5(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num);
        if(eqx(3))
            out_phi = eqx(1);
        elseif(eqx(4))
            out_phi = eqx(2);
        end
    elseif (2*pi/3 < phi_num && phi_num <= pi)
        eqx = trafo.DinY.fphi.f6(Ldab_num, n_num, Ld1_num, Ld2_num, Lm_num, Po_num, fs_num, Vi_num, Vo_num);
        if(eqx(3))
            out_phi = eqx(1);
        elseif(eqx(4))
            out_phi = eqx(2);
        end
    end
end