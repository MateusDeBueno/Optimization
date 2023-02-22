function [M,M_harmonic] = f_YY(Vp,L,n,d,fs,phi,k)
    if phi<pi/3
        M = f_YY_menor60(Vp,L,n,d,fs,phi,k);
        M_harmonic = f_YY_harmonic_menor60(Vp,L,n,d,fs,phi,k);
    else
        M = f_YY_maior60(Vp,L,n,d,fs,phi,k);
        M_harmonic = f_YY_harmonic_maior60(Vp,L,n,d,fs,phi,k);
    end
end