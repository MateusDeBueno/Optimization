clc; clear; close all;


% http://www.femm.info/list/msg01303.html

% format long
format short eng
f = 100e3; %frequency

addpath('C:\femm42\mfiles');
addpath('utils');


% https://www.mathworks.com/help/sps/ref/inverseparktransform.html


% %situacao 1
% d_q_0 = 10*[1; 0; 0];
% theta = 0*pi/3;
% fig_name = 'sit1.bmp';

% % % situacao 2
% d_q_0 = 10*[1; 0; 0];
% theta = 2*pi/3;
% fig_name = 'sit2.bmp';
% 
% %situacao 3
% d_q_0 = 10*[1; 0; 0];
% theta = 4*pi/3;
% fig_name = 'sit3.bmp';
% 
% %situacao 4
% d_q_0 = 10*[1; 0; 1];
% theta = 2*pi/3;
% fig_name = 'sit4.bmp';

%situacao 5
d_q_0 = 10*[0; 0; 1];
theta = 2*pi/3;
fig_name = 'sit5.png';

T = [sin(theta)         cos(theta)          sqrt(0.5);
    sin(theta-2*pi/3)   cos(theta-2*pi/3)    sqrt(0.5);
    sin(theta+2*pi/3)    cos(theta+2*pi/3)    sqrt(0.5)]*sqrt(2/3);

abc = T * d_q_0;



% % add circuit, series 
Np = 5; %n of turn
Ns = 5; %n of turn
Ia = abc(1)
Ib = abc(2)
Ic = abc(3)

% Ia = 10;
% Ib = 5;
% Ic = -6;


% 
% mi_addcircprop('Coil',Ia,1) ;
% mi_addcircprop('Coil',Ib,1) ;
% mi_addcircprop('Coil',Ic,1) ;




openfemm;
newdocument(0); %magnetic type problem

%gap
g = 0.1;

%carreter espessura
s = 1.25;


G_carretel = 29.2-1.6*2;
B = 29.5;

layer = 2;
N = 16*2; %numero de voltas

d_litz = 1.45;
d_w = d_litz;

s_turn = d_w*1.0095;
s_layer = d_w*0.875;

%core geometry
% G = 29.6;
% G_carretel = 29.2-1.6*2;
% E = 12.2;
% C = 42.4;
% B = 29.5;
% A = 42;
% D = 15;

G = 21.9*2;
G_carretel = G_carretel-1.6*2;
E = 22;
C = 33.2*2;
B = 48;
A = 70.5;
D = 32;



all_center = []; %salvar todos os centros de interesse, para depois integrar as regioes

mi_probdef(f,'millimeters','planar',1e-8,D);

%getting materials
mi_getmaterial('Air');
mi_getmaterial('1mm');
mi_getmaterial('0.5mm');

mi_addcircprop('Coil_a',Ia,1);
mi_addcircprop('Coil_b',Ib,1);
mi_addcircprop('Coil_c',Ic,1);

%geometrical parameters (mm)
x2 = E/2;
x3 = B/2;
x4 = A/2;
y1 = g/2;
y2 = (g+G)/2;
y3 = (g+C)/2;

%getting materials
mi_getmaterial('Air');
mi_getmaterial('Copper');


mi_addcircprop('Icoil',Ia,1) ;

%ferrite
u0 = 2100;
mu_x = u0;
mu_y = u0;
mi_addmaterial('Ferrite', mu_x, mu_y, 0);


%% CORE
horizontal_offset = -(x4+g/2);
core_draw(x2,x3,x4,y1,y2,y3,horizontal_offset)
core_draw(x2,x3,x4,y1,y2,y3,-horizontal_offset)

%air
P_air = [2*(x2+x4),0];
mi_addblocklabel(P_air);
mi_selectlabel(P_air);
mi_setblockprop('Air',0,0,'',0,0,0);
mi_clearselected();



%coil
p_coil_right = [x4-x3+s+g/2, y2-s
                x4-x3+s+g/2, -y2+s
                x4-x3+2*s+g/2, -y2+s
                x4-x3+2*s+g/2, y2-s];
DrawClosedPolygon(p_coil_right);



[p_coil_right2] = horizontal_move(p_coil_right,2*s);
DrawClosedPolygon(p_coil_right2);

p_coil_right3 = [x3+1.5*g-s, y2-s
                x3+1.5*g-s, -y2+s
                x3+1.5*g-2*s, -y2+s
                x3+1.5*g-2*s, y2-s];

DrawClosedPolygon(p_coil_right3);

[p_coil_right4] = horizontal_move(p_coil_right3,-2*s);
DrawClosedPolygon(p_coil_right4);

[p_coil_right5] = horizontal_move(p_coil_right,x3-x2+x2*2);
DrawClosedPolygon(p_coil_right5);

[p_coil_right6] = horizontal_move(p_coil_right5,2*s);
DrawClosedPolygon(p_coil_right6);

[p_coil_left] = horizontal_flip(p_coil_right);
DrawClosedPolygon(p_coil_left);

[p_coil_left2] = horizontal_flip(p_coil_right2);
DrawClosedPolygon(p_coil_left2);

[p_coil_left3] = horizontal_flip(p_coil_right3);
DrawClosedPolygon(p_coil_left3);

[p_coil_left4] = horizontal_flip(p_coil_right4);
DrawClosedPolygon(p_coil_left4);

[p_coil_left5] = horizontal_flip(p_coil_right5);
DrawClosedPolygon(p_coil_left5);

[p_coil_left6] = horizontal_flip(p_coil_right6);
DrawClosedPolygon(p_coil_left6);

%% coil fase b
P_coil = sum(p_coil_right)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_b',0,0,Np);
mi_clearselected()

P_coil = sum(p_coil_left)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_b',0,0,-Np);
mi_clearselected()

P_coil = sum(p_coil_right2)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_b',0,0,Ns);
mi_clearselected()

P_coil = sum(p_coil_left2)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_b',0,0,-Ns);
mi_clearselected()


%% coil fase c
P_coil = sum(p_coil_right5)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_c',0,0,Np);
mi_clearselected()

P_coil = sum(p_coil_right6)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_c',0,0,Ns);
mi_clearselected()

P_coil = sum(p_coil_right3)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_c',0,0,-Np);
mi_clearselected()

P_coil = sum(p_coil_right4)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_c',0,0,-Ns);
mi_clearselected()


%% coil fase a
P_coil = sum(p_coil_left4)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_a',0,0,Ns);
mi_clearselected()

P_coil = sum(p_coil_left5)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_a',0,0,-Np);
mi_clearselected()

P_coil = sum(p_coil_left3)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_a',0,0,Np);
mi_clearselected()

P_coil = sum(p_coil_left6)/4;
mi_addblocklabel(P_coil)
mi_selectlabel(P_coil);
mi_setblockprop('0.5mm',0,0,'Coil_a',0,0,-Ns);
mi_clearselected()

%%

name_fem = 'Automatic.fem';
mi_makeABC();
mi_zoomnatural();

mi_saveas(name_fem);

mi_analyze();
mi_loadsolution();

mo_showdensityplot(0,1,0.4,0,'mag')
mo_zoom(-(x4*2+s),-y3-s,+(x4*2+s),+y3+s)
mo_hidecontourplot();

mo_savebitmap(fig_name)
 



