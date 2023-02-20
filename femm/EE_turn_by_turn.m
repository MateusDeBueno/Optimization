clc; clear; close all;


% http://www.femm.info/list/msg01303.html
% calculando resistencia a partir da potencia

% format long
format short eng
% f = 100e3; %frequency
f = 100e3; %frequency
I = 1; %current

addpath('C:\femm42\mfiles');
addpath('utils');




openfemm;
newdocument(0); %magnetic type problem

%gap
g = 0.635;

%carreter espessura
s = 2;


G_carretel = 29.2-1.6*2;
B = 29.5;

layer = 2;
N = 16*2; %numero de voltas
% d_w = 1.15; %diametro de um fio awg 27
% d_w = 1.45; %diametro de um fio awg 15

d_litz = 1.45;
d_w = d_litz;

s_turn = d_w*1.0095;
s_layer = d_w*0.875;

%core geometry
G = 29.6;
G_carretel = 29.2-1.6*2;
E = 12.2;
C = 42.4;
B = 29.5;
A = 42;
D = 15;

all_center = []; %salvar todos os centros de interesse, para depois integrar as regioes

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
mi_getmaterial('Copper');

% add circuit, series 
Ia = 1;

mi_addcircprop('Icoil',Ia,1) ;

%ferrite
u0 = 2100;
mu_x = u0;
mu_y = u0;
mi_addmaterial('Ferrite', mu_x, mu_y, 0);


%% CORE
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
all_center = [all_center; P_core_u];
mi_clearselected()

p_core_d = vertical_flip(p_core_u);
DrawClosedPolygon(p_core_d);
P_core_d = [0,-y2];
mi_addblocklabel(P_core_d);
mi_selectlabel(P_core_d);
mi_setblockprop('Ferrite',0,0,'',0,0,0);
all_center = [all_center; P_core_d];
mi_clearselected()

%% COIL

% sm = G_carretel/N;
% 
% for k=1:N
%     x(k) = k*sm;
% end
% x = x - sum(x)/N;
% 
% 




position = layer_position(N,layer);
position = horizontal_multiply(position,s_layer);
position = vertical_multiply(position,s_turn);
position = horizontal_move(position,x2+s);

wires_r = position;


for k=1:length(wires_r)
    point = wires_r(k,:);
    DrawClosedCirc(point,d_w)
    mi_addblocklabel(point);
    mi_selectlabel(point);
    mi_setblockprop('Copper',0,0,'Icoil',0,0,1);
    all_center = [all_center; point];
    mi_clearselected()
end

wires_l = horizontal_flip(wires_r);
for k=1:length(wires_l)
    point = wires_l(k,:);
    DrawClosedCirc(point,d_w)
    mi_addblocklabel(point);
    mi_selectlabel(point);
    mi_setblockprop('Copper',0,0,'Icoil',0,0,-1);
    all_center = [all_center; point];
    mi_clearselected()
end


%% AIR
P_air = [x2+x4,0];
mi_addblocklabel(P_air);
mi_selectlabel(P_air);
mi_setblockprop('Air',0,0,'',0,0,0);
all_center = [all_center; P_air];
mi_clearselected();

%% CALCULATION
name_fem = 'Automatic.fem';
mi_makeABC();
mi_zoomnatural();

mi_saveas(name_fem);
mi_analyze();
mi_loadsolution();
% mo_hidecontourplot();
% 
% 
% 
% 
Coil_props = mo_getcircuitproperties('Icoil'); %currente, voltage, flux

ILa = Coil_props(1);
VLa = Coil_props(2);
Fluxa = Coil_props(3);

k = 1e3;
T = 1/f;
omega = 2*pi*f;
time = 0:T/k:T;

VLa_time = real(VLa)*cos(omega*time) - imag(VLa)*sin(omega*time);
ILa_time = real(ILa)*cos(omega*time) - imag(ILa)*sin(omega*time);
Fluxa_time = real(Fluxa)*cos(omega*time) - imag(Fluxa)*sin(omega*time);

figure
hold on
plot(time,VLa_time,'LineWidth',3)
plot(time,ILa_time*max(VLa_time),'LineWidth',3)
plot(time,VLa_time.*ILa_time,'LineWidth',3)
% plot(time,abs(VLa_time.*ILa_time),'LineWidth',3)
hold off
grid on
legend('VLa', strcat('ILa x',num2str(max(VLa_time)) ), 'PLa')

L2 = Fluxa/Ia
R = VLa/ILa
R_abs = sum(abs(VLa_time.*ILa_time))/T




