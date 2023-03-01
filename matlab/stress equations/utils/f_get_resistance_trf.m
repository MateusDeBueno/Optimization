function [Rtot] = f_get_resistance_trf(f,Np,Ns,Sp,Ss,d_l)
    
    cu = 1.72e-8; % [ohm*m]
    MLT = 230*1e-3; % [m]
    strand = 180;
    
    n = Ns/Np;
    
    Fr = f_Fr_litz2(f,1);
    
    R_dc_p = 4*cu*Np*MLT/(pi*d_l^2)/(strand*Sp);
    R_dc_s = 4*cu*Ns*MLT/(pi*d_l^2)/(strand*Ss);
    
    Rtot = Fr*(R_dc_p+R_dc_s/(n*n));

end