function L = f_get_indutance(core_geometry,N,gap,ur)

%     G_carretel = core_geometry.G_carretel;
    G = core_geometry.G*1e-3;
    E = core_geometry.E*1e-3;
    C = core_geometry.C*1e-3;
    B = core_geometry.B*1e-3;
    A = core_geometry.A*1e-3;
    D = core_geometry.D*1e-3;
    
    u0 = 4*pi*1e-7;
    % https://ieeexplore.ieee.org/document/5944575

    d = gap;
    h = G/2; %sempre

    % CENTER
    w_y = D;
    Rbasic_x_center = f_Rbasic(w_y,d,h,u0);
    w_x = E;
    Rbasic_y_center = f_Rbasic(w_x,d,h,u0);
    Rm_center = f_Rm_air_gap(d,w_x,w_y,Rbasic_x_center,Rbasic_y_center,u0);

    % SIDE
    Rbasic_x_side = f_Rbasic(w_y,d,h,u0);
    w_x = (A-B)/2;
    Rbasic_y_side = f_Rbasic(w_x,d,h,u0);
    Rm_side = f_Rm_air_gap(d,w_x,w_y,Rbasic_x_side,Rbasic_y_side,u0);


    R1_li = G/2;
    R1_Ai = D*E/2;
    [R2_li,R2_Ai] = f_corners_relutance(E/2,(C/2-G/2),h);
    R3_li = (B-E)/2;
    R3_Ai = D*(C-G)/2;
    [R4_li,R4_Ai] = f_corners_relutance((A/2-B/2),(C/2-G/2),h);
    R5_li = G/2;
    R5_Ai = D*(A-B)/2;


    R1 = R1_li/(u0*ur*R1_Ai);
    R2 = R2_li/(u0*ur*R2_Ai);
    R3 = R3_li/(u0*ur*R3_Ai);
    R4 = R4_li/(u0*ur*R4_Ai);
    R5 = R5_li/(u0*ur*R5_Ai);

    Rt_side = R3+R4+R5+Rm_side+R5+R4+R3;
    Rt_center = R2+R1+Rm_center+R1+R2;
    Rt = Rt_center + Rt_side*0.5;

    L = N*N/Rt;
end

