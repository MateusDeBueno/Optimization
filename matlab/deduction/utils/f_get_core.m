function core_geometry = f_get_core(code)
    %core geometry
    if code==1
        core_name = 422115;
        G = 14.8*2;
        G_carretel = G-1.6*2;
        E = 12.2;
        C = 21.2*2;
        B = 29.5;
        A = 42;
        D = 15.5;
        Ve = 17.6e-6;
        Ae = 181e-6; 
        Ac = E*1e-3*D*1e-3;
        Rth = 19;
    elseif code == 2
        core_name = 422120;
        G = 14.8*2;
        G_carretel = G-1.6*2;
        E = 12.2;
        C = 21.2*2;
        B = 29.5;
        A = 42;
        D = 20;
        Ve = 23.3e-6;
        Ae = 240e-6; 
        Ac = E*1e-3*D*1e-3;
        Rth = 15;
    elseif code == 3
        core_name = 552821;
        G = 18.5*2;
        G_carretel = G-1.6*2;
        E = 17.2;
        C = 27.8*2;
        B = 37.5;
        A = 55;
        D = 21;
        Ve = 42.5e-6;
        Ae = 354e-6; 
        Ac = E*1e-3*D*1e-3;
        Rth = 11;
    elseif code == 4
        core_name = 552825;
        G = 18.5*2;
        G_carretel = G-1.6*2;
        E = 17.2;
        C = 27.8*2;
        B = 37.5;
        A = 55;
        D = 25;
        Ve = 52.1e-6;
        Ae = 420e-6; 
        Ac = E*1e-3*D*1e-3;
        Rth = 8;
    end
    core_geometry.core_name = core_name;
    core_geometry.G_carretel = G_carretel;core_geometry.G = G;core_geometry.E = E;core_geometry.C = C;core_geometry.B = B;
    core_geometry.A = A;core_geometry.D = D; core_geometry.Ve = Ve; core_geometry.Ae = Ae; core_geometry.Ac = Ac;
    core_geometry.Rth = Rth;
end

