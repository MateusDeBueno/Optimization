clc; clear; close all;

% format long
format short eng
f = 1000e3; %frequency
I = 0.1; %current

addpath('C:\femm42\mfiles');
addpath('utils');

openfemm;
newdocument(0); %magnetic type problem

%gap
g = 0.65;

%carreter espessura
s = 2;
space = 2;

%core geometry
G = 29.6;
G_carretel = 29.2-1.6*2;
E = 12.2;
C = 42.4;
B = 29.5;
A = 42;
D = 20;

mi_probdef(f,'millimeters','planar',1e-8,D);

%geometrical parameters (mm)
x2 = E/2;
x3 = B/2;
x4 = A/2;
y1 = g/2;
y2 = (g+G)/2;
y3 = (g+C)/2;

%getting materials
mi_getmaterial('Air');
mi_getmaterial('1mm');
mi_getmaterial('0.5mm');

% add circuit, series 
Ia = 1;
N = 57; %n of turn
mi_addcircprop('Coil',Ia,1) ;

%ferrite
u0 = 2100;
mu_x = u0;
mu_y = u0;
mi_addmaterial('Ferrite', mu_x, mu_y, 0);

p_core_u_right = [x2,y1
                 x2,y2
                 x3,y2
                 x3,y1
                 x4,y1
                 x4,y3];
p_core_u = [p_core_u_right; flip(horizontal_flip(p_core_u_right))];
DrawClosedPolygon(p_core_u);
P_core_u = [0,y2];
mi_addblocklabel(P_core_u);
mi_selectlabel(P_core_u);
mi_setblockprop('Ferrite',0,0,'',0,0,0);
mi_clearselected()

p_core_d = vertical_flip(p_core_u);
DrawClosedPolygon(p_core_d);
P_core_d = [0,-y2];
mi_addblocklabel(P_core_d);
mi_selectlabel(P_core_d);
mi_setblockprop('Ferrite',0,0,'',0,0,0);
mi_clearselected()

%coild
p_coil_right = [x2+s, y2-s
                x3-space, y2-s
                x3-space, -y2+s
                x2+s, -y2+s];
DrawClosedPolygon(p_coil_right);
P_coil_r = sum(p_coil_right)/4;
mi_addblocklabel(P_coil_r)
mi_selectlabel(P_coil_r);
mi_setblockprop('0.5mm',0,0,'Coil',0,0,N);
mi_clearselected()
p_coil_left = horizontal_flip(p_coil_right);
DrawClosedPolygon(p_coil_left);
P_coil_l = sum(p_coil_left)/4;
mi_addblocklabel(P_coil_l);
mi_selectlabel(P_coil_l);
mi_setblockprop('0.5mm',0,0,'Coil',0,0,-N);
mi_clearselected();

%air
P_air = [x2+x4,0];
mi_addblocklabel(P_air);
mi_selectlabel(P_air);
mi_setblockprop('Air',0,0,'',0,0,0);
mi_clearselected();

name_fem = 'Automatic.fem';
mi_makeABC();
mi_zoomnatural();

mi_saveas(name_fem);
mi_analyze();
mi_loadsolution();
mo_hidecontourplot();
mo_showdensityplot(1,0,1.25e-3,0,'bmag');

% mo_selectblock(P_air)
% mo_selectblock(P_core_u)
% mo_selectblock(P_core_d)
% mo_selectblock(P_coil_r)
% mo_selectblock(P_coil_l)
% 
% w = mo_blockintegral(2);
% L1 = 2*w/(Ia*Ia)



Coil_props = mo_getcircuitproperties('Coil'); %currente, voltage, flux
L2 = Coil_props(3)/Ia
R100k = Coil_props(2)/Coil_props(1)
