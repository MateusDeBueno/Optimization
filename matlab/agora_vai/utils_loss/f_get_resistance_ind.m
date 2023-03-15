function Rind = f_get_resistance_ind(h,n,l)
    %h: harmonic
    %n: porosity
    %l: dab struct

    cu = 1.72e-8; % [ohm*m]
    
    dl = l.L.dl;
    strand = l.L.strand;
    Fr = f_Fr_litz(h,l.L.Nl,dl,strand,n,l);
    R_s_dc = 4*cu*l.L.N*l.L.MLT/(pi*l.L.dl^2);
    R_w_dc = R_s_dc/l.L.strand; %R_DC

    Rind = Fr*R_w_dc;
end