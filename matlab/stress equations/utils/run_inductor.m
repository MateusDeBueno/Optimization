function [Rtot,Ltot,n_layers,position] = run_inductor(core_geometry,wire_geometry,g,f)

    G_carretel = core_geometry.G_carretel;
    G = core_geometry.G;
    E = core_geometry.E;
    C = core_geometry.C;
    B = core_geometry.B;
    A = core_geometry.A;
    D = core_geometry.D;
    s = core_geometry.s;
    
    d_litz = wire_geometry.d_litz;
    strands = wire_geometry.strands;
    sbw = wire_geometry.sbw;
    N = wire_geometry.N;
    
    Ia = 2;

    x2 = E/2;
    x3 = B/2;
    x4 = A/2;
    y1 = g/2;
    y2 = (g+G)/2;
    y3 = (g+C)/2;

    [position_litz,position,group,max_spires_by_layer] = layer_position_litz2(d_litz,strands,N,sbw,G_carretel,x2,s);


    % FEMM INITIALIZATION

    openfemm(1);
    newdocument(0); %magnetic type problem
    hidepointprops;
    mi_probdef(f,'millimeters','planar',1e-8,D); %config unit

    %getting materials
    mi_getmaterial('Air');
    mi_getmaterial('Copper');

    %ferrite
    u0 = 2100;
    mu_x = u0;
    mu_y = u0;
    mi_addmaterial('Ferrite', mu_x, mu_y, 0);

    %create circuits
    for k=1:N
        mi_addcircprop(strcat('Coil_',num2str(k),'r'),Ia,0)
        mi_addcircprop(strcat('Coil_',num2str(k),'l'),-Ia,0)
    end

    % CORE
    core_draw(x2,x3,x4,y1,y2,y3);

    % COIL
    %plot all wires, in the same string is parallell
    put_wires(position_litz,d_litz,1,group,'r')
    put_wires(horizontal_flip(position_litz),d_litz,-1,group,'l')

    % AIR

    air_draw(x2,x4)

    % CALCULATION
    name_fem = 'Automatic.fem';
    mi_makeABC();
    mi_zoomnatural();

    mi_saveas(name_fem);
    mi_analyze();
    mi_loadsolution();


    % GET VALUES

    values_right = zeros(3,N);
    values_left = zeros(3,N);
    amp_right = zeros(1,N);
    amp_left = zeros(1,N);
    voltage_right = zeros(1,N);
    voltage_left = zeros(1,N);
    resistance_right = zeros(1,N);
    resistance_left = zeros(1,N);
    flux_right = zeros(1,N);
    flux_left = zeros(1,N);

    for k=1:N
        values_right(:,k) = mo_getcircuitproperties(strcat('Coil_',num2str(k),'r'));
        values_left(:,k) = mo_getcircuitproperties(strcat('Coil_',num2str(k),'l'));

        amp_right(k) = values_right(1,k);
        amp_left(k) = values_left(1,k);

        voltage_right(k) = values_right(2,k);
        voltage_left(k) = values_left(2,k);

        flux_right(k) = values_right(3,k);
        flux_left(k) = values_left(3,k);

        resistance_right(k) = voltage_right(k)/amp_right(k);
        resistance_left(k) = voltage_left(k)/amp_left(k);

        flux_right(k) = flux_right(k)/amp_right(k);
        flux_left(k) = flux_left(k)/amp_left(k);    
    end

    Rtot = sum(real(resistance_right))+sum(real(resistance_left));
    Ltot = sum(real(flux_right))+sum(real(flux_left));
    n_layers = ceil(N/max_spires_by_layer);
end

