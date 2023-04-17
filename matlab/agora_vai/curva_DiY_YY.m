clear; close all; clc;

%% iDY

phi_exp{1,1} = [-15 -19];
io_exp{1,1} = [6.41 5.169]; %phi 10 20 30 40 50
ef_exp{1,1} = [88 85.81]/100;

phi_exp{1,2} = [-8 -14 -20 -24];
io_exp{1,2} = [8.252 6.437 4.607 3.381]; %phi 10 20 30 40 50
ef_exp{1,2} = [93.874 92.621 90.409 87.648]/100;

phi_exp{1,3} = [-27 -22 -18 -15 -12 -8 -5 -2];
io_exp{1,3} = [2.478 3.972 5.217 6.208 7.07 8.3 9.285 10.112]; %phi 10 20 30 40 50
ef_exp{1,3} = [88.1 92.04 93.583 94.348 94.811 95.336 95.696 95.745]/100;

phi_exp{1,4} = [-31 -28 -25 -21 -18 -15 -11 -8];
io_exp{1,4} = [1.085 1.958 2.929 4.14 5 5.963 7.17 8]; %phi 10 20 30 40 50
ef_exp{1,4} = [82.7 89.72 92.99 94.39 95.4 95.99 96.38 96.58]/100;

phi_exp{1,5} = [-28 -25 -22 -19 -15 -12 -9 -5 -2 0];
io_exp{1,5} = [1.141 2.119 2.971 3.943 5.142 5.983 6.941 8.135 8.956 9.410]; %phi 10 20 30 40 50
ef_exp{1,5} = [90.541 94.509 95.845 96.645 97.126 97.315 97.44 97.545 97.599 97.552]/100;

phi_exp{1,6} = [-22 -20 -18 -15 -12 -8 -5 -2 1 4];
io_exp{1,6} = [1.9 2.5 3.15 4.123 5.295 6.15 7.168 7.999 9.086 10]; %phi 10 20 30 40 50
ef_exp{1,6} = [93.166 94.624 95.55 96.33 96.8 96.96 97.23 97.33 97.356 97.28]/100;

%% YY

phi_exp{2,1} = [nan];
io_exp{2,1} = [nan]; %phi 10 20 30 40 50
ef_exp{2,1} = [nan]/100;

phi_exp{2,2} = [10 20 30 40 50 58];
io_exp{2,2} = [2.34 3.824 5.34 6.77 8.04 8.918]; %phi 10 20 30 40 50
ef_exp{2,2} = [95.2 97.3 97.4 97.3 97.1 96.8]/100;

phi_exp{2,3} = [nan];
io_exp{2,3} = [nan]; %phi 10 20 30 40 50
ef_exp{2,3} = [nan]/100;

phi_exp{2,4} = [10 20 30 40 50 60 70 75];
io_exp{2,4} = [1 2.88 4.56 6.13 7.647 8.88 9.77 10]; %phi 10 20 30 40 50
ef_exp{2,4} = [80.9 92.2 94.5 95.82 96.38 96.134 95.78 95.45]/100;

phi_exp{2,5} = [nan];
io_exp{2,5} = [nan]; %phi 10 20 30 40 50
ef_exp{2,5} = [nan]/100;

phi_exp{2,6} = [10 20 30 40 50 60 70 79];
io_exp{2,6} = [.813 2.68 4.366 5.887 7.23 8.404 9.539 10]; %phi 10 20 30 40 50
ef_exp{2,6} = [70.04 88.36 92.077 93.65 94.27 94.595 94.6 93.99]/100;


%%

ocultar_titulo = 0;

% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue

efc_min_graph = 0.9;

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

l.L.a = 1.285676564102240;
l.L.b = 2.349349709620920;
l.L.kc = 11.187048860336922;
% l.L.a = 1.394;
% l.L.b = 2.248;
% l.L.kc = 2.448;
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

% l.tr.a = 1.585;
% l.tr.b = 2.756;
% l.tr.kc = 0.402;
l.tr.a = 1.448345520625372;
l.tr.b = 2.729356903304897;
l.tr.kc = 2.925783179535760;
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

l.sC.Rac100 = 4.3e-3 + 2e-3;
l.C.Rb_ac10k = 4.2e-3 + 2e-3;

l.pr.dt = 0e-9;
l.pr.Vi = 400;
l.pr.d = 1;
l.pr.fs = 100e3;
l.pr.Ldab = 61.5e-6;
l.pr.Ld1 = 1.4e-6;
l.pr.n = l.tr.Ns/l.tr.Np;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 1500e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);

%% calcula iDY

pr_ang = 30; %precision
ang_min(1) = -29.5*pi/180;
ang_max(1) = 60*pi/180;
vec_ph(1,:) = ang_min(1):(ang_max(1)-ang_min(1))/pr_ang:ang_max(1);

ang_min(2) = .5*pi/180;
ang_max(2) = 90*pi/180;
vec_ph(2,:) = ang_min(2):(ang_max(2)-ang_min(2))/pr_ang:ang_max(2);

pr_d = 30; %precision
d_max = 1.2;
d_min = 0.25;
vec_d = d_min:0.125/2:d_max;
vec_d = [vec_d;vec_d];
n_points = length(vec_d)*length(vec_ph(1,:));

trafo_usados = {'DiY', 'YY'};%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

for nf = 1:2
    k = 0;
    trafoo = trafo_usados(nf);
    for kk=1:length(vec_ph(nf,:))
        for ii = 1:length(vec_d)
            k = k + 1;
    
            phii(nf,k) = vec_ph(nf,kk);
            dd(nf,k) = vec_d(nf,ii);
    
            l.pr.phi = phii(nf,k);
            l.pr.d = dd(nf,k);
            
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
            
            [efc(nf,k),ptot(nf,k),Pm(nf,k),cSw_p(nf,k),cSw_s(nf,k),sSw_p(nf,k),sSw_s(nf,k),Pil_cu(nf,k),PiL_cu(nf,k),PiLd_cu(nf,k),Psc(nf,k),P_core_L(nf,k),P_core_tr(nf,k)] = f_loss(trafoo, out, l);
%             output{nf,k} = out;
        end
    end
end

%% eficiencia em funcao da tensao de saida e do angulo
contourLevels = [0.93 0.95 0.97];
colors = f_create_cmap(length(contourLevels)+1, color1, color2);

for nf = 1:2
    y = dd(nf,:);
    x = rad2deg(phii(nf,:));
    z = efc(nf,:);
    
    n = 200;
    %Create regular grid across data space
    [X,Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
    Z = griddata(x,y,z,X,Y);
    
    fig = figure;
    hold on
    set(fig,'defaultLegendAutoUpdate','off');
    
    [legen_name] = create_legend_contourf(contourLevels, colors);
    hc = contourfcmap(X,Y*l.pr.Vi,Z,contourLevels, colors(2:end-1,:), ...
         'lo', colors(1,:), ...
         'hi', colors(end,:), ...
         'method', 'calccontour');
    hc.h.LineStyle = 'none';
    hold off
    legend(legen_name,'Location','best','FontSize', 16,'Interpreter','latex');
    
    ylim([min(vec_d(:))*l.pr.Vi max(vec_d(:))*l.pr.Vi])
    
    grid on
    grid minor
    set(gca, 'FontSize', 20)
    ylabel('$V_o\,$[V]')
    xlabel('$\phi\,[^{\circ}]$')
    if (ocultar_titulo == 0)
        title(trafo_usados(nf))
    end
    file_name = append('figure\comparacoes_perdas\mapa_eficiencia_Vo_phi_',trafo_usados{nf},'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');    
end

%% eficia em funcao da tensao de saida e potencia de saida

contourLevels = [0.93 0.95 0.97];
colors = f_create_cmap(length(contourLevels)+1, color1, color2);

for nf = 1:2
    y = dd(nf,:);
    x = Pm(nf,:);
    z = efc(nf,:);
    
    n = 200;
    %Create regular grid across data space
    [X,Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
    Z = griddata(x,y,z,X,Y);


    % Find unique values of dd
    unique_dd = unique(dd(nf,:));
    
    % Loop through unique values of dd and find the maximum value of Pm for each
    max_Pm = zeros(size(unique_dd));
    for i = 1:length(unique_dd)
        idx = (dd(nf,:) == unique_dd(i));
        max_Pm(i) = max(Pm(idx));
    end
    
    fig = figure;
    hold on
    set(fig,'defaultLegendAutoUpdate','off');
    
    [legen_name] = create_legend_contourf(contourLevels, colors);
    hc = contourfcmap(X*1e-3,Y*l.pr.Vi,Z,contourLevels, colors(2:end-1,:), ...
         'lo', colors(1,:), ...
         'hi', colors(end,:), ...
         'method', 'calccontour');
    hc.h.LineStyle = 'none';
    
    hold off
    legend(legen_name,'Location','best','FontSize', 16,'Interpreter','latex');
    
    ylim([min(vec_d(:))*l.pr.Vi max(vec_d(:))*l.pr.Vi])
    
    grid on
    grid minor
    set(gca, 'FontSize', 20)
    ylabel('$V_o\,$[V]')
    xlabel('$P_o\,$[kW]')

    if (ocultar_titulo == 0)
        title(trafo_usados(nf))
    end

    file_name = append('figure\comparacoes_perdas\mapa_eficiencia_Vo_Po_',trafo_usados{nf},'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
end

% % % % %% analise com tensoes fixas part1
% % % % 
% % % % d_values = [150/400 200/400 250/400];
% % % % 
% % % % % Initialize arrays to store values
% % % % phi_v = cell(2, numel(d_values));
% % % % Pm_v = cell(2, numel(d_values));
% % % % efc_v = cell(2, numel(d_values));
% % % % 
% % % % colors = f_create_cmap(length(d_values), color1, color2);
% % % % for nf = 1:2
% % % % 
% % % %     clear legen_name;
% % % %     % Loop over d values
% % % %     for i = 1:numel(d_values)
% % % %         % Get logical index for current d value
% % % %         idx = round(dd(nf,:), 3) == d_values(i);
% % % %         sum(idx)
% % % %         % Extract values for current d value
% % % %         phi_v{nf,i} = phii(nf,idx);
% % % %         Pm_v{nf,i} = Pm(nf,idx);
% % % %         efc_v{nf,i} = efc(nf,idx);
% % % %     
% % % %         legen_name(i) = append('$V_o = $',string(d_values(i)*l.pr.Vi),'$\,$[V]');
% % % %     end
% % % % end
% % % % 
% % % % for nf=1:2
% % % %     figure
% % % %     hold on
% % % %     for i = 1:numel(d_values)
% % % %         plot(Pm_v{nf,i}./(l.pr.Vi*d_values(i)), efc_v{nf,i}, 'LineWidth', 1.5, 'Color', colors(i,:))
% % % %     end
% % % %     hold off
% % % %     ylim([efc_min_graph 1])
% % % %     xlim([0 10])
% % % %     legend(legen_name, 'Location','best')
% % % %     set(gca, 'FontSize', 20)
% % % %     grid on
% % % %     grid minor
% % % %     ylabel('$\eta$')
% % % %     xlabel('$i_o\,$[A]')
% % % %     if (ocultar_titulo == 0)
% % % %         title(trafo_usados(nf))
% % % %     end    
% % % % 
% % % %     figure
% % % %     hold on
% % % %     for i = 1:numel(d_values)
% % % %         plot(Pm_v{nf,i}*1e-3, efc_v{nf,i}, 'LineWidth', 1.5, 'Color', colors(i,:))
% % % %     end
% % % %     hold off
% % % %     ylim([efc_min_graph 1])
% % % %     xlim([0 4])
% % % %     legend(legen_name, 'Location','best')
% % % %     set(gca, 'FontSize', 20)
% % % %     grid on
% % % %     grid minor
% % % %     ylabel('$\eta$')
% % % %     xlabel('$P_o\,$[kW]') 
% % % % 
% % % %     if (ocultar_titulo == 0)
% % % %         title(trafo_usados(nf))
% % % %     end
% % % % end
% % % % 
% % % % 
% % % % %% analise com tensoes fixas part2
% % % % 
% % % % d_values = [300/400 350/400 400/400];
% % % % 
% % % % % Initialize arrays to store values
% % % % phi_v = cell(2, numel(d_values));
% % % % Pm_v = cell(2, numel(d_values));
% % % % efc_v = cell(2, numel(d_values));
% % % % 
% % % % colors = f_create_cmap(length(d_values), color1, color2);
% % % % for nf = 1:2
% % % % 
% % % %     clear legen_name;
% % % %     % Loop over d values
% % % %     for i = 1:numel(d_values)
% % % %         % Get logical index for current d value
% % % %         idx = round(dd(nf,:), 3) == d_values(i);
% % % %         sum(idx)
% % % %         % Extract values for current d value
% % % %         phi_v{nf,i} = phii(nf,idx);
% % % %         Pm_v{nf,i} = Pm(nf,idx);
% % % %         efc_v{nf,i} = efc(nf,idx);
% % % %     
% % % %         legen_name(i) = append('$V_o = $',string(d_values(i)*l.pr.Vi),'$\,$[V]');
% % % %     end
% % % % end
% % % % 
% % % % for nf=1:2
% % % %     figure
% % % %     hold on
% % % %     for i = 1:numel(d_values)
% % % %         plot(Pm_v{nf,i}./(l.pr.Vi*d_values(i)), efc_v{nf,i}, 'LineWidth', 1.5, 'Color', colors(i,:))
% % % %     end
% % % %     hold off
% % % %     ylim([efc_min_graph 1])
% % % %     xlim([0 10])
% % % %     legend(legen_name, 'Location','best')
% % % %     set(gca, 'FontSize', 20)
% % % %     grid on
% % % %     grid minor
% % % %     ylabel('$\eta$')
% % % %     xlabel('$i_o\,$[A]')
% % % %     if (ocultar_titulo == 0)
% % % %         title(trafo_usados(nf))
% % % %     end
% % % %     
% % % %     figure
% % % %     hold on
% % % %     for i = 1:numel(d_values)
% % % %         plot(Pm_v{nf,i}*1e-3, efc_v{nf,i}, 'LineWidth', 1.5, 'Color', colors(i,:))
% % % %     end
% % % %     hold off
% % % %     ylim([efc_min_graph 1])
% % % %     xlim([0 4])
% % % %     legend(legen_name, 'Location','best')
% % % %     set(gca, 'FontSize', 20)
% % % %     grid on
% % % %     grid minor
% % % %     ylabel('$\eta$')
% % % %     xlabel('$P_o\,$[kW]') 
% % % %     if (ocultar_titulo == 0)
% % % %         title(trafo_usados(nf))
% % % %     end    
% % % % end

%% analisando com tensoes fixas, par a par

d_values = [150/400 200/400 250/400 300/400 350/400 400/400];

% Initialize arrays to store values
phi_v = cell(2, numel(d_values));
Pm_v = cell(2, numel(d_values));
efc_v = cell(2, numel(d_values));

colors = f_create_cmap(length(d_values), color1, color2);
for nf = 1:2
    clear legen_name;
    % Loop over d values
    for i = 1:numel(d_values)
        % Get logical index for current d value
        idx = round(dd(nf,:), 3) == d_values(i);
        
        % Extract values for current d value
        phi_v{nf,i} = phii(nf,idx);
        Pm_v{nf,i} = Pm(nf,idx);
        efc_v{nf,i} = efc(nf,idx);
    end
end

for d_target = d_values
    colors = f_create_cmap(2, color1, color2);
    idxx = d_target == d_values;
    position = find(idxx == 1);
    figure
    hold on
    for nf = 1:2
        plot(Pm_v{nf,position}./(l.pr.Vi*d_target), efc_v{nf,position}, 'LineWidth', 1.5, 'Color', colors(nf,:))
    end
    for nf = 1:2
        scatter(io_exp{nf,position},ef_exp{nf,position},'filled','square','MarkerEdgeColor',colors(nf,:),'MarkerFaceColor',colors(nf,:))
    end

    hold off
    legend(trafo_usados, 'Location','best')
    ylim([efc_min_graph 1])
    xlim([0 10])
    set(gca, 'FontSize', 20)
    grid on
    grid minor
    ylabel('$\eta$')
    xlabel('$i_o\,$[A]')
    if (ocultar_titulo == 0)
        title(append('$V_o = $',string(d_target*l.pr.Vi),'$\,$[V]'))
    end
    file_name = append('figure\comparacoes_perdas\perdas_io_Vo_',string(d_target*l.pr.Vi),'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
end


%% todas as perdas

d_values = [150/400 200/400 250/400 300/400 350/400 400/400];

% Initialize arrays to store values
phi_v = cell(2, numel(d_values));
Pm_v = cell(2, numel(d_values));
efc_v = cell(2, numel(d_values));

cSw_p_v = cell(2, numel(d_values));
cSw_s_v = cell(2, numel(d_values));
sSw_p_v = cell(2, numel(d_values));
sSw_s_v = cell(2, numel(d_values));

Pil_cu_v = cell(2, numel(d_values));
PiL_cu_v = cell(2, numel(d_values));
PiLd_cu_v = cell(2, numel(d_values));

P_core_L_v = cell(2, numel(d_values));
P_core_tr_v = cell(2, numel(d_values));




colors = f_create_cmap(length(d_values), color1, color2);
for nf = 1:2
    clear legen_name;
    % Loop over d values
    for i = 1:numel(d_values)
        % Get logical index for current d value
        idx = round(dd(nf,:), 3) == d_values(i);
        sum(idx)
        % Extract values for current d value
        phi_v{nf,i} = phii(nf,idx);
        Pm_v{nf,i} = Pm(nf,idx);
        efc_v{nf,i} = efc(nf,idx);
        
        cSw_p_v{nf,i} = cSw_p(nf,idx); %perdas conducao primario
        cSw_s_v{nf,i} = cSw_s(nf,idx); %perdas conducao secundario
        sSw_p_v{nf,i} = sSw_p(nf,idx); %perdas comutacao primario
        sSw_s_v{nf,i} = sSw_s(nf,idx); %perdas comutacao secundario
        
        Pil_cu_v{nf,i} = Pil_cu(nf,idx); %perdas bobina primario
        PiL_cu_v{nf,i} = PiL_cu(nf,idx); %perdas bobina secundario
        PiLd_cu_v{nf,i} = PiLd_cu(nf,idx); %perdas bobina indutor
        
        P_core_L_v{nf,i} = P_core_L(nf,idx); %perdas nucleo indutor
        P_core_tr_v{nf,i} = P_core_tr(nf,idx); %perdas nucleo trafo
        
        ptot_v{nf,i} = ptot(nf,idx); %perdas totais
        % legen_name(i) = append('$V_o = $',string(d_values(i)*l.pr.Vi),'$\,$[V]');
    end
end

legen_name = {'iDY: $Sw_p$', 'iDY: $Sw_s$','YY: $Sw_p$', 'YY: $Sw_s$'};
for d_target = d_values
    colors = f_create_cmap(2, color1, color2);
    idxx = d_target == d_values;
    position = find(idxx == 1);
    fig = figure;
    hold on
    set(fig,'defaultLegendAutoUpdate','off');
    
    for nf = 1:2
        plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (cSw_p_v{nf,position}+sSw_p_v{nf,position}),':o', 'LineWidth', 1.5, 'Color', colors(nf,:))
        plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (cSw_s_v{nf,position}+sSw_s_v{nf,position}),'--.', 'LineWidth', 1.5, 'Color', colors(nf,:))
%         plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (PiLd_cu_v{nf,position}+P_core_L_v{nf,position}),'--o', 'LineWidth', 1.5, 'Color', colors(nf,:))
%         plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (Pil_cu_v{nf,position}+PiL_cu_v{nf,position}+P_core_tr_v{nf,position}),'--*', 'LineWidth', 1.5, 'Color', colors(nf,:))
%         plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (ptot_v{nf,position}),'-', 'LineWidth', 1.5, 'Color', colors(nf,:))
    end
    hold off
    legend(legen_name, 'Location','best','NumColumns',2)
%     legend(legen_name, 'Location','northwest','NumColumns',2)
    xlim([0 10])
    set(gca, 'FontSize', 20)
    grid on
    grid minor
    ylabel('$P\,$[W]')
    xlabel('$i_o\,$[A]')
    if (ocultar_titulo == 0)
        title(append('$V_o = $',string(d_target*l.pr.Vi),'$\,$[V]'))
    end
    file_name = append('figure\comparacoes_perdas\perdas_chaves_Vo_',string(d_target*l.pr.Vi),'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
end


legen_name = {'iDY:$L_{dab}$', 'iDY:Trans','YY:$L_{dab}$', 'YY:Trans'};
for d_target = d_values
    colors = f_create_cmap(2, color1, color2);
    idxx = d_target == d_values;
    position = find(idxx == 1);
    fig = figure;
    hold on
    set(fig,'defaultLegendAutoUpdate','off');

    for nf = 1:2
%         plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (cSw_p_v{nf,position}+sSw_p_v{nf,position}),':o', 'LineWidth', 1.5, 'Color', colors(nf,:))
%         plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (cSw_s_v{nf,position}+sSw_s_v{nf,position}),'--.', 'LineWidth', 1.5, 'Color', colors(nf,:))
        plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (PiLd_cu_v{nf,position}+P_core_L_v{nf,position}),'--o', 'LineWidth', 1.5, 'Color', colors(nf,:))
        plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (Pil_cu_v{nf,position}+PiL_cu_v{nf,position}+P_core_tr_v{nf,position}),'--*', 'LineWidth', 1.5, 'Color', colors(nf,:))
%         plot(Pm_v{nf,position}./(l.pr.Vi*d_target), (ptot_v{nf,position}),'-', 'LineWidth', 1.5, 'Color', colors(nf,:))
    end
    hold off
    legend(legen_name, 'Location','best','NumColumns',2)
%     legend(legen_name, 'Location','northwest','NumColumns',2)
    xlim([0 10])
    set(gca, 'FontSize', 20)
    grid on
    grid minor
    ylabel('$P\,$[W]')
    xlabel('$i_o\,$[A]')
    if (ocultar_titulo == 0)
        title(append('$V_o = $',string(d_target*l.pr.Vi),'$\,$[V]'))
    end
    file_name = append('figure\comparacoes_perdas\perdas_magneticos_Vo_',string(d_target*l.pr.Vi),'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
end
