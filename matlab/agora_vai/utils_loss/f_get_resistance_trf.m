function [R_ac_p, R_ac_s] = f_get_resistance_trf(h,n,l)
    
    cu = 1.72e-8; % [ohm*m]
    
    dl = l.tr.dl;
    strand = l.tr.strand;
    Fr = f_Fr_litz(h,l.L.Nl,dl,strand,n,l);
    
    R_dc_p = 4*cu*l.tr.Np*l.tr.MLT/(pi*dl^2)/(strand*l.tr.Sp);
    R_dc_s = 4*cu*l.tr.Ns*l.tr.MLT/(pi*dl^2)/(strand*l.tr.Ss);
    
    R_ac_p = Fr*R_dc_p;
    R_ac_s = Fr*R_dc_s;
end