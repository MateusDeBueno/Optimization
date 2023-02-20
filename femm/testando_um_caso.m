clc; clear; close all;

format short eng

addpath('C:\femm42\mfiles');
addpath('utils');
addpath('data');

awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32

%parametros gerais
L_goal = 65e-6; %indutancia desejada
Imax = 12; %max current, used for minimum turns
Bsat = 0.25; %usado para numero minimo de espiras
f = 100e3; %frequency
A_min = 1; %mm quadrados de cobre, usado para calcular numero d strands
A_max = 4;
Pcore_max = 15;
strands_max_total = 30;

w = 2*pi*f;
Ts = 1/f;
d = 1;
phi = 70*pi/180;
n = 15/27;
t_db = 277e-9;
T__DB_N = t_db/Ts;
Vp = 400;
IN = Vp/(w*L_goal);



%parametros do ferrite
alpha = 1.394;
beta = 2.248;
kc = 2.448;
ki = 0.183005238;

alpha = 1.64;
beta = 2.737;
kc = 0.248;
ki = 0.011;

ur = 2200;


sbw = 1.02;
i = 1;

N = 14;



core_geometry = f_get_core(4);
    
core_name = core_geometry.core_name;
G_carretel = core_geometry.G_carretel;
G = core_geometry.G;
E = core_geometry.E;
C = core_geometry.C;
B = core_geometry.B;
A = core_geometry.A;
D = core_geometry.D;
Ve = core_geometry.Ve;
Ae = core_geometry.Ae;
Ac = core_geometry.Ac;
Rth = core_geometry.Rth;
s = 2; %carretel espessura
core_geometry.s = s;

% Nmin = ceil(Imax*L_goal/(Bsat*Ac))
Imax = N*core_geometry.Ac*Bsat/L_goal;

strands = 320;
awg_chosen = 38;
d_litz = awg_wires(awg_chosen);

d_target = f_get_d_target(strands,d_litz,sbw);


wire_geometry.N = N;
wire_geometry.d_litz = d_litz;
wire_geometry.strands = strands;
wire_geometry.sbw = sbw;
wire_geometry.awg_chosen = awg_chosen;


g = f_get_gap(core_geometry,N,ur,L_goal)
output = f_trafo_YY(d,phi,n,T__DB_N,alpha,Vp,Ae,N,f);

Irsm = output(1,10)*IN;
I_L_max = output(1,end-1)*IN;
integrais_L = output(1,end);

deltaB = 2*L_goal*I_L_max/(N*Ae);
watts_per_m3 = ki/(Ts/2)*abs(deltaB)^(beta-alpha)*integrais_L;
Pcore = Ve*watts_per_m3

%GET RESISTANCE AND INDUCTANCE
[Rtot,Ltot,vector_n_layers,vector_position] = run_inductor(core_geometry,wire_geometry,g,f);


%WINDING LOSSES
Pwinding = Rtot*Irsm^2;

%TOTAL LOSSES
Ptot = Pwinding + Pcore;

%TEMPERATURE
Temp = Rth*Pcore


