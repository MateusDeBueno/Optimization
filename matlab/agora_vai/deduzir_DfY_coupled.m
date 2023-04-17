clear; close all; clc;

% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

addpath('utils')
addpath('utils_transf')
addpath('utils_loss')
addpath('data');
addpath('contourfcmap-pkg-master/FEX-function_handle')
addpath('contourfcmap-pkg-master/arclength')
addpath('contourfcmap-pkg-master/contourcs')
addpath('contourfcmap-pkg-master/contourfcmap')
addpath('contourfcmap-pkg-master/distance2curve')
addpath('contourfcmap-pkg-master/fillnan')
addpath('contourfcmap-pkg-master/interparc')
addpath('contourfcmap-pkg-master/minmax')
addpath('contourfcmap-pkg-master/multiplepolyint')
addpath('contourfcmap-pkg-master/parsepv')
addpath('contourfcmap-pkg-master/pcolorbar')
color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue

syms Ld1 Ld2 n Lm Po t L1 L2 Ldab M fs Vi d dt real positive
syms phi real

Vo = Vi*d;
Ts = 1/fs;

phi_num = [deg2rad(-10),deg2rad(50)];  %[MUDAR]
[x0s.f1, ts.f1, idab.f1, hb.f1, HB.f1, Ip.f1, Is.f1, iME.f1, idrm.f1, ilrm.f1, iLrm.f1, iSwPrm.f1, iSwSrm.f1] = simplify_DfY(phi_num(1));
[x0s.f2, ts.f2, idab.f2, hb.f2, HB.f2, Ip.f2, Is.f2, iME.f2, idrm.f2, ilrm.f2, iLrm.f2, iSwPrm.f2, iSwSrm.f2] = simplify_DfY(phi_num(2));

trafo = 'DfY';

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61.5e-6;
Ld1_num = 2e-6;
n_num = 1;
Ld2_num = 2e-6;
Lm_num = [.5e-3, 1.5e-3 10e-3];
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

%cria intervalo de angulo
intervalo.f1 = [deg2rad(-30) 0];
intervalo.f2 = [0 deg2rad(60)];

%% passo para plot
pr = 200; %precision
vec_ph.f2 = min(intervalo.f2):(max(intervalo.f2)-min(intervalo.f2))/pr:max(intervalo.f2);
vec_ph.f1 = min(intervalo.f1):(max(intervalo.f1)-min(intervalo.f1))/pr:max(intervalo.f1);


fidab.f1 = matlabFunction(idab.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fidab.f2 = matlabFunction(idab.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fx0s.f1 = matlabFunction(x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fx0s.f2 = matlabFunction(x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fts.f1 = matlabFunction(ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fts.f2 = matlabFunction(ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%equacao_limite
Ip_eq.f1 = Ip.f1 == 0;
Is_eq.f1 = Is.f1 == 0;
pot_eq.f1 = iME.f1*Vo == Po;

fpot_eq.f1 = matlabFunction(solve(pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIp_eq.f1 = matlabFunction(solve(Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs_eq.f1 = matlabFunction(solve(Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

Ip_eq.f2 = Ip.f2 == 0;
Is_eq.f2 = Is.f2 == 0;
pot_eq.f2 = iME.f2*Vo == Po;

fpot_eq.f2 = matlabFunction(solve(pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIp_eq.f2 = matlabFunction(solve(Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs_eq.f2 = matlabFunction(solve(Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%correntes eficazes
fidrm.f1 = matlabFunction(idrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
filrm.f1 = matlabFunction(ilrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiLrm.f1 = matlabFunction(iLrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwPrm.f1 = matlabFunction(iSwPrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwSrm.f1 = matlabFunction(iSwSrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

fidrm.f2 = matlabFunction(idrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
filrm.f2 = matlabFunction(ilrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiLrm.f2 = matlabFunction(iLrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwPrm.f2 = matlabFunction(iSwPrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fiSwSrm.f2 = matlabFunction(iSwSrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

fIp.f1 = matlabFunction(Ip.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs.f1 = matlabFunction(Is.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIp.f2 = matlabFunction(Ip.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
fIs.f2 = matlabFunction(Is.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

mag=1;

fIp.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1))
fIs.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1))
fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1))

fIp.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(2))
fIs.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(2))
fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(2))

%% corrente nos estados
% 
% yy = fx0s.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1));
% xx = fts.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1));
% 
% %corrente no primario do trafo
% figure
% cmap = f_create_cmap(3, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% hold on
% 
% plot(xx*1e6,yy(1,:),'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(xx*1e6,yy(2,:),'Color',jetcustom(2,:),'LineWidth',1.5)
% plot(xx*1e6,yy(3,:),'Color',jetcustom(3,:),'LineWidth',1.5)
% 
% hold off
% grid on
% grid minor
% set(gca, 'FontSize', 20)
% xlabel('$t[\mu$s]')
% ylabel('$i\,$[A]')
% legend({'$i_{a}$','$i_{b}$','$i_{c}$'},'Location','southeast','FontSize', 16)
% file_name = append('figure\finalCap2\primary_current_',trafo,'.pdf');
% exportgraphics(gca,file_name,'ContentType','vector');
% 
% 
% %corrente no secundario do trafo
% figure
% cmap = f_create_cmap(3, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% hold on
% 
% plot(xx*1e6,yy(4,:),'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(xx*1e6,yy(5,:),'Color',jetcustom(2,:),'LineWidth',1.5)
% plot(xx*1e6,yy(6,:),'Color',jetcustom(3,:),'LineWidth',1.5)
% 
% hold off
% grid on
% grid minor
% set(gca, 'FontSize', 20)
% xlabel('$t[\mu$s]')
% ylabel('$i\,$[A]')
% legend({'$i_{A}$','$i_{B}$','$i_{C}$'},'Location','southeast','FontSize', 16)
% file_name = append('figure\finalCap2\secondary_current_',trafo,'.pdf');
% exportgraphics(gca,file_name,'ContentType','vector');
% 
% % corrente no Ldab
% yy_s = fx0s.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1));
% xx = fts.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(1));
% il = yy_s(1:3,:);
% yy = il - [il(3,:); il(1:2,:)]; %equacao especifica para esse trafo
% 
% figure
% cmap = f_create_cmap(3, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% hold on
% plot(xx*1e6,yy(1,:),'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(xx*1e6,yy(2,:),'Color',jetcustom(2,:),'LineWidth',1.5)
% plot(xx*1e6,yy(3,:),'Color',jetcustom(3,:),'LineWidth',1.5)
% hold off
% grid on
% grid minor
% set(gca, 'FontSize', 20)
% xlabel('$t[\mu$s]')
% ylabel('$i\,$[A]')
% legend({'$i_{Ldab_a}$','$i_{Ldab_b}$','$i_{Ldab_c}$'},'Location','southeast','FontSize', 16)

%%










%% corrente de comutacao

contourLevels = [-3 0 3];
colors = f_create_cmap(length(contourLevels)+1, color1, color2);

Xf1 = vec_ph.f1;
Xf2 = vec_ph.f2;

vec_d = 0:0.1:3;
for ii = 1:length(vec_d)
    for k = 1:length(vec_ph.f1)

        dd = vec_d(ii);
        phii = vec_ph.f1(k);

        Zpp.f1(k,ii) = fIp.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,dd,fs_num,phii);
        Zss.f1(k,ii) = fIs.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,dd,fs_num,phii);
    end
end
for ii = 1:length(vec_d)
    for k = 1:length(vec_ph.f2)

        dd = vec_d(ii);
        phii = vec_ph.f2(k);

        Zpp.f2(k,ii) = fIp.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,dd,fs_num,phii);
        Zss.f2(k,ii) = fIs.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,dd,fs_num,phii);
    end
end

%analisando aqui
fIs.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,1,fs_num,0)

y = vec_d;
x = Xf1*180/pi;
z = Zpp.f1';
n = 200;
[X1,Y1] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z1 = griddata(x,y,z,X1,Y1);

y = vec_d;
x = Xf2*180/pi;
z = Zpp.f2';
n = 200;
[X2,Y2] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z2 = griddata(x,y,z,X2,Y2);

fig = figure;
set(fig,'defaultLegendAutoUpdate','off');
hold on

[legen_name] = create_legend_contourf(contourLevels, colors);

leg = legend(legen_name,'Location','southeast','FontSize', 16,'Interpreter','latex');
title(leg,'$I_{sw-p}\,$[A]')

hold on
hc = contourfcmap(X1,Y1,Z1,contourLevels, colors(2:end-1,:), ...
     'lo', colors(1,:), ...
     'hi', colors(end,:), ...
     'method', 'calccontour');
hc.h.LineStyle = 'none';

hc = contourfcmap(X2,Y2,Z2,contourLevels, colors(2:end-1,:), ...
     'lo', colors(1,:), ...
     'hi', colors(end,:), ...
     'method', 'calccontour');
hc.h.LineStyle = 'none';
hold off

grid on
grid minor
set(gca, 'FontSize', 20)
ylabel('d')
xlabel('$\phi\,[^{\circ}]$')

file_name = append('figure\finalCap2\Ip_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

y = vec_d;
x = Xf1*180/pi;
z = Zss.f1';
n = 200;
[X1,Y1] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z1 = griddata(x,y,z,X1,Y1);

y = vec_d;
x = Xf2*180/pi;
z = Zss.f2';
n = 200;
[X2,Y2] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z2 = griddata(x,y,z,X2,Y2);

fig = figure;
set(fig,'defaultLegendAutoUpdate','off');
hold on

[legen_name] = create_legend_contourf(contourLevels, colors);
leg = legend(legen_name,'Location','best','FontSize', 16,'Interpreter','latex');
title(leg,'$I_{sw-s}\,$[A]')

hold on
hc = contourfcmap(X1,Y1,Z1,contourLevels, colors(2:end-1,:), ...
     'lo', colors(1,:), ...
     'hi', colors(end,:), ...
     'method', 'calccontour');
hc.h.LineStyle = 'none';

hc = contourfcmap(X2,Y2,Z2,contourLevels, colors(2:end-1,:), ...
     'lo', colors(1,:), ...
     'hi', colors(end,:), ...
     'method', 'calccontour');
hc.h.LineStyle = 'none';
hold off

grid on
grid minor
set(gca, 'FontSize', 20)
ylabel('d')
xlabel('$\phi\,[^{\circ}]$')

file_name = append('figure\finalCap2\Is_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% região de ZVS

cor_aqui = color2;
c1 = f_create_cmap(3, cor_aqui, [1 1 1]);
color_white = [c1(2,:)];

mag = 1;
figure
hold on

h = fill([1 1], [1 1],color_white,'Edgecolor', 'none');
h.EdgeColor = color2;
h.LineWidth = 1.5;

lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);

fill([rad2deg(vec_ph.f1) fliplr(rad2deg(vec_ph.f1))], [lim_p.f1 fliplr(lim_s.f1)],color_white, 'EdgeColor', 'none');
fill([rad2deg(vec_ph.f2) fliplr(rad2deg(vec_ph.f2))], [lim_p.f2 fliplr(lim_s.f2)],color_white, 'EdgeColor', 'none');


plot(rad2deg(vec_ph.f1),lim_p.f1,'Color',cor_aqui,'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_p.f2,'Color',cor_aqui,'LineWidth',1.5)
plot(rad2deg(vec_ph.f1),lim_s.f1,'Color',cor_aqui,'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_s.f2,'Color',cor_aqui,'LineWidth',1.5)

text(-20,2.5,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(-20,0.45,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(30,1.5,'SS in both','Interpreter', 'Latex','FontSize', 14) 

hold off
ylim([0 3])
grid on
grid minor
xlabel('$\phi\,[^{\circ}]$')
ylabel('$d$')
set(gca, 'FontSize', 20)
xlim([min(intervalo.f1)*180/pi max(intervalo.f2)*180/pi])
file_name = append('figure\finalCap2\ZVS_phi_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

% plot d por Po
figure
hold on

lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
fill([lim_pot.f1 fliplr(lim_pot.f1)], [lim_p.f1 fliplr(lim_s.f1)],color_white, 'EdgeColor', 'none');
lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
fill([lim_pot.f2 fliplr(lim_pot.f2)], [lim_p.f2 fliplr(lim_s.f2)],color_white, 'EdgeColor', 'none');
plot(lim_pot.f1,lim_p.f1,'Color',cor_aqui,'LineWidth',1.5)
plot(lim_pot.f2,lim_p.f2,'Color',cor_aqui,'LineWidth',1.5)
plot(lim_pot.f1,lim_s.f1,'Color',cor_aqui,'LineWidth',1.5)
plot(lim_pot.f2,lim_s.f2,'Color',cor_aqui,'LineWidth',1.5)

text(250,2.5,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(250,0.45,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(1000,1.5,'SS in both','Interpreter', 'Latex','FontSize', 14) 

hold off
ylim([0 3])
xlim([0 max(lim_pot.f2)])
grid on
grid minor
xlabel('$P_o$[W]')
ylabel('$d$')
set(gca, 'FontSize', 20)
file_name = append('figure\finalCap2\ZVS_Po_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');







%% plot Po por phi
mag = 1;
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
hold on

lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),lim_pot.f1,'Color','k','LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_pot.f2,'Color','k','LineWidth',1.5)

hold off
ylim([0 1500])
xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
grid on
grid minor
xlabel('$\phi\,[^{\circ}]$')
ylabel('$P_o\,$[W]')
set(gca, 'FontSize', 20)

file_name = append('figure\finalCap2\gain_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');





%% potencia de saida

contourLevels = [250 750 1000 1250]/1000;
colors = f_create_cmap(length(contourLevels)+1, color1, color2);

Xf1 = vec_ph.f1;
Xf2 = vec_ph.f2;

vec_d = 0:0.1:3;
for ii = 1:length(vec_d)
    for k = 1:length(vec_ph.f1)

        dd = vec_d(ii);
        phii = vec_ph.f1(k);
        Poo.f1(k,ii) = fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,dd,fs_num,phii);
    end
end
for ii = 1:length(vec_d)
    for k = 1:length(vec_ph.f2)

        dd = vec_d(ii);
        phii = vec_ph.f2(k);

        Poo.f2(k,ii) = fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,dd,fs_num,phii);
    end
end


y = vec_d;
x = Xf1*180/pi;
z = Poo.f1';
n = 200;
[X1,Y1] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z1 = griddata(x,y,z,X1,Y1);

y = vec_d;
x = Xf2*180/pi;
z = Poo.f2';
n = 200;
[X2,Y2] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z2 = griddata(x,y,z,X2,Y2);
hc = contourfcmap(X1,Y1,Z1/1000,contourLevels, colors(2:end-1,:), ...
     'lo', colors(1,:), ...
     'hi', colors(end,:), ...
     'method', 'calccontour');
hc.h.LineStyle = 'none';
file_name = append('figure\finalCap2\Po_map_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');


%% indutancia mutua

% phi_chosen = max(intervalo.f1);
phi_chosen = deg2rad(15);
for mag=1:3
    equat.f2 = fx0s.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_chosen);
    equat_ts.f2 = fts.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_chosen);
    
    figure
    cmap = f_create_cmap(2, color2, color1);
    colormap(cmap)
    jetcustom = cmap;
    hold on
    plot(equat_ts.f2*1e6,equat.f2(1,:),'Color',jetcustom(1,:),'LineWidth',1.5)
    plot(equat_ts.f2*1e6,-equat.f2(4,:), '-.','Color',jetcustom(2,:),'LineWidth',1.5)
    hold off
    xlim([0 max(equat_ts.f2)*1e6])
    grid on
    grid minor
    legend({'$i_{tra}$','$-i_{trA}$'},'Location','best','FontSize', 14)
    set(gca, 'FontSize', 20)
    xlabel('$t[\mu$s]')
    ylabel('$i\,$[A]')
    ylim([-3 3])
    
    file_name = append('figure\finalCap2\iaiA_Lm_mag',num2str(mag),'_',trafo,'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');    
end
%%
LLm_num = min(Lm_num):(max(Lm_num)-min(Lm_num))/100:max(Lm_num);
MM_num = LLm_num*n_num;
LL1_num = Ld1_num + LLm_num;
LL2_num = Ld2_num + n_num*n_num*LLm_num;

for k = 1:length(vec_ph.f2)
    for ii = 1:length(LLm_num)

        phi_ii = vec_ph.f2(k);
        Lm_ii = LLm_num(ii);
        MM_num = Lm_ii*n_num;
        LL1_num = Ld1_num + Lm_ii;
        LL2_num = Ld2_num + n_num*n_num*Lm_ii;     

        v_ilrm.f2(k,ii) = filrm.f2(LL1_num,LL2_num,Ldab_num,MM_num,Vi_num,d_num,fs_num,phi_ii);
        v_iLrm.f2(k,ii) = fiLrm.f2(LL1_num,LL2_num,Ldab_num,MM_num,Vi_num,d_num,fs_num,phi_ii);
    end
end

for k = 1:length(vec_ph.f1)
    for ii = 1:length(LLm_num)

        phi_ii = vec_ph.f1(k);
        Lm_ii = LLm_num(ii);
        MM_num = Lm_ii*n_num;
        LL1_num = Ld1_num + Lm_ii;
        LL2_num = Ld2_num + n_num*n_num*Lm_ii;     

        v_ilrm.f1(k,ii) = filrm.f1(LL1_num,LL2_num,Ldab_num,MM_num,Vi_num,d_num,fs_num,phi_ii);
        v_iLrm.f1(k,ii) = fiLrm.f1(LL1_num,LL2_num,Ldab_num,MM_num,Vi_num,d_num,fs_num,phi_ii);
    end
end

x1 = LLm_num;
y1 = rad2deg(vec_ph.f1);
z1 = v_iLrm.f1./v_ilrm.f1;

x2 = LLm_num;
y2 = rad2deg(vec_ph.f2);
z2 = v_iLrm.f2./v_ilrm.f2;

%% com fill

contourLevels = [0.94 0.99 1 1.01];
colors = f_create_cmap(length(contourLevels)+1, color1, color2);

fig = figure;
set(fig,'defaultLegendAutoUpdate','off');
hold on

[legen_name] = create_legend_contourf(contourLevels, colors);

leg = legend(legen_name,'Location','best','FontSize', 14,'Interpreter','latex');

hold on

x = x1*1e3;
y = y1;
z = z1;
n = 200;
[X1,Y1] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z1 = griddata(x,y,z,X1,Y1);

hc = contourfcmap(X1,Y1,Z1,contourLevels, colors(2:end-1,:), ...
     'lo', colors(1,:), ...
     'hi', colors(end,:), ...
     'method', 'calccontour');
hc.h.LineStyle = 'none';

x = x2*1e3;
y = y2;
z = z2;
n = 200;
[X2,Y2] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
Z2 = griddata(x,y,z,X2,Y2);

hc = contourfcmap(X2,Y2,Z2,contourLevels, colors(2:end-1,:), ...
     'lo', colors(1,:), ...
     'hi', colors(end,:), ...
     'method', 'calccontour');
hc.h.LineStyle = 'none';

hold off

grid on
grid minor
set(gca, 'FontSize', 20)
xlabel('$L_m\,$[mH]')
ylabel('$\phi\,[^{\circ}]$')

file_name = append('figure\finalCap2\LM_ratio_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');


%%
% Create the contour plot
figure
hold on
[C,h] = contour(x1*1e3, y1, z1,contourLevels, 'ShowText', 'on','Color',color2);
h.LineWidth = 1.5;
clabel(C,h,'FontSize',15,'Color',[0 0 0],'LineWidth',2)
[C,h] = contour(x2*1e3, y2, z2,contourLevels, 'ShowText', 'on','Color',color2);
h.LineWidth = 1.5;
clabel(C,h,'FontSize',15,'Color',[0 0 0],'LineWidth',2)
hold off
grid on
grid minor
xlabel('$L_m\,$[mH]')
ylabel('$\phi\,[^{\circ}]$')
set(gca, 'FontSize', 20)
file_name = append('figure\finalCap2\LM_ratio_contour_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');





%% plot d por phi
figure
hold on

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
    lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
    lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
    lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
    fill([rad2deg(vec_ph.f1) fliplr(rad2deg(vec_ph.f1))], [lim_p.f1 fliplr(lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    fill([rad2deg(vec_ph.f2) fliplr(rad2deg(vec_ph.f2))], [lim_p.f2 fliplr(lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(rad2deg(vec_ph.f1),lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(rad2deg(vec_ph.f2),lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(rad2deg(vec_ph.f1),lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(rad2deg(vec_ph.f2),lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end

text(-20,2.5,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(-20,0.45,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(30,2.5,'SS in both','Interpreter', 'Latex','FontSize', 14) 

hold off
ylim([0 3])
grid on
grid minor

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

xlabel('$\phi\,[^{\circ}]$')
ylabel('$d$')
set(gca, 'FontSize', 20)
xlim([min(intervalo.f1)*180/pi max(intervalo.f2)*180/pi])
% f_save_figure(append('figure\',string(trafo),'_d_phi.pdf'))

file_name = append('figure\finalCap2\ZVS_phi_Lm_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%% plot d por Po


cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

c1 = f_create_cmap(3, cmap(1,:), [1 1 1]);
c2 = f_create_cmap(3, cmap(2,:), [1 1 1]);
c3 = f_create_cmap(3, cmap(3,:), [1 1 1]);
color_white = [c1(2,:);c2(2,:);c3(2,:)];

figure
hold on


for ii=1:3
    h = fill([1 1], [1 1],color_white(ii,:),'Edgecolor', 'none');
    h.EdgeColor = jetcustom(ii,:);
    h.LineWidth = 1.5;
end

for mag=1:3
    lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
    lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
    lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f1);
    fill([lim_pot.f1 fliplr(lim_pot.f1)], [lim_p.f1 fliplr(lim_s.f1)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
    lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
    lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,vec_ph.f2);
    fill([lim_pot.f2 fliplr(lim_pot.f2)], [lim_p.f2 fliplr(lim_s.f2)],color_white(mag,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
    plot(lim_pot.f1,lim_p.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(lim_pot.f2,lim_p.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(lim_pot.f1,lim_s.f1,'Color',jetcustom(mag,:),'LineWidth',1.5)
    plot(lim_pot.f2,lim_s.f2,'Color',jetcustom(mag,:),'LineWidth',1.5)
end

text(300,2.5,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(800,0.35,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(1000,1.55,'SS in both','Interpreter', 'Latex','FontSize', 14) 

hold off
ylim([0 3])
xlim([0 max(lim_pot.f2)])
grid on
grid minor

for iN = 1:length(Lm_num)
    legendCell{iN} = append('$L_m =$ ',num2str(Lm_num(iN)*1000),'$\,$mH');
end
legend(legendCell,'Location','best','FontSize', 14)

xlabel('$P_o$[W]')
ylabel('$d$')
set(gca, 'FontSize', 20)

file_name = append('figure\finalCap2\ZVS_Lm_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');
