function [ph] = YY_get_angle(trafo,precision,Ptarget)

    [k] = get_k(vP,Ptarget);
    ph = vec_phi(k);
end