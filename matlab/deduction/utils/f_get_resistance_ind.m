function Rind = f_get_resistance_ind(f,N,d_l)
    cu = 1.72e-8; % [ohm*m]
    
    MLT = 140*1e-3; % [m]
    strand = 180;
    
    Fr = f_Fr_litz2(f,1);
    
    R_s_dc = 4*cu*N*MLT/(pi*d_l^2);
    R_w_dc = R_s_dc/strand; %R_DC
    
    
    Rind = Fr*R_w_dc;
end