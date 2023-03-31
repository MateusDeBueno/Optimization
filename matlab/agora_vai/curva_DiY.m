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

awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32

load('YY.mat');
load('YD.mat');
load('DfD.mat');
load('DiD.mat');
load('DiY.mat');
load('DfY.mat');

% Combine the structs into a single struct
l.eq = struct('YY', YY, 'YD', YD, 'DfD', DfD, ...
              'DiD', DiD, 'DiY', DiY, 'DfY', DfY);

l.sw = load('f_fitted_off.mat');
l.sw = load('f_fitted_on.mat');
l.sw.Ronp = 112e-3;
l.sw.Rons = 112e-3;


l.L.a = 1.394;
l.L.b = 2.248;
l.L.kc = 2.448;
l.L.int_ki = integral(@(theta) abs(cos(theta)).^l.L.a,0,2*pi);
l.L.ki = l.L.kc/(2^(l.L.b-l.L.a)*(2*pi)^(l.L.a-1)*l.L.int_ki);
l.L.Ac = 421.3e-6;
l.L.Ve = 52.1e-6;
l.L.N = 14;
l.L.strand = 180;
l.L.awg = 38;
l.L.dl = awg_wires(l.L.awg)*1e-3; % [m]
l.L.MLT = 140*1e-3; % [m]
l.L.Nl = 1;

l.tr.a = 1.585;
l.tr.b = 2.756;
l.tr.kc = 0.402;
l.tr.int_ki = integral(@(theta) abs(cos(theta)).^l.tr.a,0,2*pi);
l.tr.ki = l.tr.kc/(2^(l.tr.b-l.tr.a)*(2*pi)^(l.tr.a-1)*l.tr.int_ki);
l.tr.Ac = 683e-6;
l.tr.Ve = 146.79e-6;
l.tr.Np = 9;
l.tr.Ns = 5;
l.tr.Sp = 2;
l.tr.Ss = 3;
l.tr.N = l.tr.Ns;
l.tr.strand = 180;
l.tr.awg = 38;
l.tr.dl = awg_wires(l.tr.awg)*1e-3; % [m]
l.tr.MLT = 230*1e-3; % [m]

l.sC.Rac100 = 4.3e-3;
l.C.Rb_ac10k = 4.2e-3;


l.pr.dt = 250e-9;
% l.pr.phi = deg2rad(67);
l.pr.Vi = 400;
l.pr.d = 1;
l.pr.fs = 100e3;
l.pr.Ldab = 61e-6;
l.pr.Ld1 = 1.4e-6;
l.pr.n = l.tr.Ns/l.tr.Np;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 700e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);

%% TRAFOOOOOOOOOOOOOO
trafoo = 'DiY';

%%

color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue

pr_ang = 200; %precision
ang_min = -29*pi/180;
ang_max = 60*pi/180;
vec_ph = ang_min:(ang_max-ang_min)/pr_ang:ang_max;
pr_d = 20; %precision
d_max = 1.2;
d_min = 0.2;
vec_d = d_min:(d_max-d_min)/pr_d:d_max;

n_points = length(vec_d)*length(vec_ph)

k = 0;
for kk=1:length(vec_ph)
    for ii = 1:length(vec_d)
        k = k + 1;

        phii(k) = vec_ph(kk);
        dd(k) = vec_d(ii);

        l.pr.phi = phii(k);
        l.pr.d = dd(k);

        if (trafoo == "YY")
            out = f_equationsYY_dt(l);
        elseif (trafoo == "YD")
            out = f_equationsYD_dt(l);
        elseif (trafoo == "DfD")
            out = f_equationsDfD_dt(l);
        elseif (trafoo == "DiD")
            out = f_equationsDiD_dt(l);
        elseif (trafoo == "DiY")
            out = f_equationsDiY_dt(l);
        elseif (trafoo == "DfY")
            out = f_equationsDfY_dt(l);
        end

        [efc(k),ptot(k),Pm(k),cSw_p(k),cSw_s(k),sSw_p(k),sSw_s(k),Pil_cu(k),PiL_cu(k),PiLd_cu(k),Psc(k),P_core_L(k),P_core_tr(k)] = f_loss(trafoo, out, l);
    end
end

%%

y = dd;
x = rad2deg(phii);
z = efc;

n = 200;
%Create regular grid across data space
[X,Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));



figure
contourLevels = [0 0.9 0.92 0.94 0.96];
colors = f_create_cmap(length(contourLevels), color1, color2);
[cCont, hCont] = contourf(X,Y*l.pr.Vi,griddata(x,y,z,X,Y),contourLevels, 'LineStyle', 'none');
cmap = interp1(contourLevels, colors, linspace(min(contourLevels), max(contourLevels), 256));
colormap(cmap);
leg = contourLegend(hCont);
leg.FontSize = 14;
grid on
grid minor
set(gca, 'FontSize', 20)
ylabel('$V_o\,$[V]')
xlabel('$\phi\,[^{\circ}]$')
file_name = append('figure\losses\Vo_phi_efc_',string(trafoo),'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');

%%

% % Find unique values of dd
% unique_dd = unique(dd);
% 
% % Loop through unique values of dd and find the maximum value of Pm for each
% max_Pm = zeros(size(unique_dd));
% for i = 1:length(unique_dd)
%     idx = (dd == unique_dd(i));
%     max_Pm(i) = max(Pm(idx));
% end
% 
% y = dd;
% x = Pm;
% z = efc;
% n = 200;
% 
% 
% %Create regular grid across data space
% [X,Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
% 
% figure
% hold on
% contourLevels = [0 0.9 0.92 0.94 0.96];
% colors = f_create_cmap(length(contourLevels), color1, color2);
% [~, hCont] = contourf(X,Y*l.pr.Vi,griddata(x,y,z,X,Y),contourLevels, 'LineStyle', 'none');
% cmap = interp1(contourLevels, colors, linspace(min(contourLevels), max(contourLevels), 256));
% colormap(cmap);
% leg = contourLegend(hCont);
% leg.FontSize = 14;
% leg.AutoUpdate = 'off';
% 
% %p max
% plot(max_Pm, unique_dd*l.pr.Vi,'LineWidth',1.5, 'Color', [0 0 0])
% xlabel('dd')
% ylabel('Max Pm')
% allChildren = get(gca, 'Children');                % list of all objects on axes
% displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
% delete(allChildren(strcmp(displayNames, 'data1')))
% hold off
% grid on
% grid minor
% ylabel('$V_o\,$[V]')
% xlabel('$P_o\,$[W]')
% set(gca, 'FontSize', 20)
% text(1.05*max_Pm(ceil(length(unique_dd)/2)),unique_dd(ceil(length(unique_dd)/2))*l.pr.Vi,'$\leftarrow\, P_{o-limit}$','Interpreter', 'Latex','FontSize', 14) 
% f_save_figure(append('figure\losses\Vo_P_efc_',string(trafoo),'.pdf'))



%% analise com tensoes fixas

phi_exp_200 = [-8 -14 -20 -24 ];
io_exp_200 = [8.252 6.437 4.607 3.381]; %phi 10 20 30 40 50
ef_exp_200 = [93.874 92.621 90.409 87.648]/100;

phi_exp_300 = [-31 -28 -25 -21 -18 -15 -11 -8];
io_exp_300 = [1.085 1.958 2.929 4.14 5 5.963 7.17 8]; %phi 10 20 30 40 50
ef_exp_300 = [82.7 89.72 92.99 94.39 95.4 95.99 96.38 96.58]/100;

phi_exp_400 = [-22 -20 -18 -15 -12 -8 -5 -2 1 4];
io_exp_400 = [1.9 2.5 3.15 4.123 5.295 6.15 7.168 7.999 9.086 10]; %phi 10 20 30 40 50
ef_exp_400 = [93.166 94.624 95.55 96.33 96.8 96.96 97.23 97.33 97.356 97.28]/100;


d_values = [0.5, 0.75, 1];

% Initialize arrays to store values
phi_v = cell(1, numel(d_values));
Pm_v = cell(1, numel(d_values));
efc_v = cell(1, numel(d_values));

clear legen_name;
% Loop over d values
for i = 1:numel(d_values)
    % Get logical index for current d value
    idx = round(dd, 2) == d_values(i);
    
    % Extract values for current d value
    phi_v{i} = phii(idx);
    Pm_v{i} = Pm(idx);
    efc_v{i} = efc(idx);

    legen_name(i) = append('$V_o = $',string(d_values(i)*l.pr.Vi),'$\,$[V]');
end


colors = f_create_cmap(3, color1, color2);

figure
hold on
for i = 1:1:numel(d_values)
    plot(Pm_v{i}./(l.pr.Vi*d_values(i)), efc_v{i}, 'LineWidth', 1.5, 'Color', colors(i,:))
end

scatter(io_exp_200,ef_exp_200,'filled','square','MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:))
scatter(io_exp_300,ef_exp_300,'filled','square','MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:))
scatter(io_exp_400,ef_exp_400,'filled','square','MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:))


grid on
ylim([0.85 1])
xlim([0 10])
legend(legen_name, 'Location','best')
set(gca, 'FontSize', 20)
grid on
grid minor
ylabel('$\eta$')
xlabel('$i_o\,$[A]')
file_name = append('figure\losses\efc_Io_',string(trafoo),'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');
%% analise das perdas

% pswitch = cSw_p + cSw_s + sSw_p + sSw_s;
% 
% y = dd;
% x = Pm;
% z = pswitch;
% 
% n = 200;
% %Create regular grid across data space
% [X,Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
% 
% 
% figure
% contourLevels = [0 20 40 60 80 100];
% colors = f_create_cmap(length(contourLevels), color1, color2);
% [~, hCont] = contourf(X,Y*l.pr.Vi,griddata(x,y,z,X,Y),contourLevels, 'LineStyle', 'none');
% cmap = interp1(contourLevels, colors, linspace(min(contourLevels), max(contourLevels), 256));
% colormap(cmap);
% leg = contourLegend(hCont);
% leg.FontSize = 14;
% grid on
% grid minor

%%

% Find unique values of dd
unique_dd = unique(dd);

% Loop through unique values of dd and find the maximum value of Pm for each
max_Pm = zeros(size(unique_dd));
for i = 1:length(unique_dd)
    idx = (dd == unique_dd(i));
    max_Pm(i) = max(Pm(idx));
end

y = dd;
x = Pm;
z = efc;

[f_fitted, goodness_off] = fit([x', y'],z','cubicinterp');

[X,Y]=meshgrid(x',y');
Z = f_fitted(X,Y);

figure
hold on
contourLevels = [0 0.9 0.92 0.94 0.96];
colors = f_create_cmap(length(contourLevels), color1, color2);
[cCont, hCont] = contourf(X,Y*l.pr.Vi,Z,contourLevels, 'LineStyle', 'none');
cmap = interp1(contourLevels, colors, linspace(min(contourLevels), max(contourLevels), 256));
colormap(cmap);
leg = contourLegend(hCont);
leg.FontSize = 14;
leg.AutoUpdate = 'off';

%p max
plot(max_Pm, unique_dd*l.pr.Vi,'LineWidth',1.5, 'Color', [0 0 0])
xlabel('dd')
ylabel('Max Pm')
allChildren = get(gca, 'Children');                % list of all objects on axes
displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
delete(allChildren(strcmp(displayNames, 'data1')))
hold off
grid on
grid minor
ylabel('$V_o\,$[V]')
xlabel('$P_o\,$[W]')
xlim([0 5000])
set(gca, 'FontSize', 20)
text(1.05*max_Pm(ceil(length(unique_dd)/2)),unique_dd(ceil(length(unique_dd)/2))*l.pr.Vi,'$\leftarrow\, P_{o-limit}$','Interpreter', 'Latex','FontSize', 14) 
file_name = append('figure\losses\Vo_P_efc_',string(trafoo),'.pdf');
exportgraphics(gca,file_name,'Resolution',500)
