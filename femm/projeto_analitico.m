clear; close all; clc




%dados do conversor
Vp = 400;
fs = 100e3;
n = 15/27;
w = fs*2*pi;
t_db = 277e-9;
% t_db = 0;
Ts = 1/fs;
T__DB_N = t_db/Ts;

%dados do ferrite
alpha = 1.22;
beta = 2.848;
ki = 4.539113451;



% http://www.femm.info/list/msg01303.html
% calculando resistencia a partir da potencia

global u0
global cu

u0 = 4*pi*1e-7;
cu = 1.72e-8;

format short eng
f = 100e3; %frequency

addpath('C:\femm42\mfiles');
addpath('utils');
addpath('data');

%max current
Imax = 10;
Bsat = 0.3;

%corrente rms na fundamental
Irms = 6;

%indutancia desejada
L_goal = 65e-6;

%core geometry
core_name = '42/21/15';
G_carretel = 29.2-1.6*2;
G = 29.6;
E = 12.2;
C = 42.4;
B = 29.5;
A = 42;
D = 15;
s = 1; %carretel espessura
MLT = 87;
Ve = 17.6e-6;
Ae = 181; %mm^2
ur = 2100;
Ac = E*1e-3*D*1e-3;

awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32


awg_chosen = 27;
d_litz = awg_wires(awg_chosen);
strands = 20;
sbw = 1.05;
p = 9.5e-6; %in m
core_geometry.G_carretel = G_carretel;core_geometry.G = G;core_geometry.E = E;core_geometry.C = C;core_geometry.B = B;core_geometry.A = A;core_geometry.D = D;core_geometry.s = s;
wire_geometry.d_litz = d_litz;wire_geometry.strands = strands;wire_geometry.sbw = sbw;




d = 1;
phi = 70;
n = 15/27;
Ipico = 10;
Vp = 400;
T__DB_N = 0.02;

w = fs*2*pi;
IN = Vp/(w*L_goal);

%parametros do ferrite
alpha = 1.22;
beta = 2.848;
ki = 4.539113451;

A_min = 1.5; %entre 1.5 e 3 mm quadrados de cobre
A_max = 3;

i = 1;



for awg_chosen=26:32 %chosen wire
   
    
    
    d_litz = awg_wires(awg_chosen);

    A_litz = pi*(d_litz/2)^2;
    
    strands_min = ceil(A_min/A_litz);
    strands_max = floor(A_max/A_litz);
%     d_target = f_get_d_target(strands,d_litz,sbw);
%     N_by_layer_max = floor(G_carretel/d_target);
%     layer_max = floor((0.5*(B-E)-s)/d_target);
% 
%     Nmin = ceil(Imax*L_goal/(Bsat*Ac));
%     Nmax = N_by_layer_max*layer_max;
% 
%     if Nmax>50
%         Nmax = 50;
%     end
    

    if (strands_max>40)
        strands_max = 40;
    end

    for strands = strands_min:strands_max
        
        d_target = f_get_d_target(strands,d_litz,sbw);
        N_by_layer_max = floor(G_carretel/d_target);
        layer_max = floor((0.5*(B-E)-s)/d_target);
        
        Nmin = ceil(Imax*L_goal/(Bsat*Ac));
        Nmax = N_by_layer_max*layer_max;
        
        if (Nmin < N_by_layer_max)
            Nmin = N_by_layer_max;
        end
        
        for N=Nmin:Nmax
            wire_geometry.N = N;

            %GAT CALCULATION
            g = f_get_gap(core_geometry,N,ur,L_goal);

            Nl = ceil(N/N_by_layer_max); %n de layers

            Acu = d_litz*strands; %in mm2

            Rdc = MLT*1e-3*N*cu/(Acu*(1e-3)^2);
            Fr100k = get_Fr(f,d_litz*1e-3,strands,Nl,p);


            output{i} = f_trafo_YY(d,phi*pi/180,n,T__DB_N,alpha,Vp,Ae*1e-9,N,f);


            ILrms = output{i}(10)*IN;
            Pcu(i) = Fr100k*Rdc*ILrms^2;
            Ipico = output{i}(22)*IN;

            integrais = output{i}(23);
            P_fe(i) = f_perdas_indutor_ferrite(L_goal,Ipico,f,alpha,beta,ki,Ve*1e-9,N,Ae*1e-9,integrais);
            i = i + 1;
        end
    end
end
figure
plot(P_fe)

figure
plot(Pcu)

figure
hold on
plot(P_fe)
plot(Pcu)
hold off
ylim([0 150])