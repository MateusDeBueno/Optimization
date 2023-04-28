clear; close all; clc;
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

l.sC.Rac100 = (4.3e-3)/2;
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
%% primeiro trafo oDY

d_vetor = [1 1 0.5 1 1 0.5];
phi_vetor = [deg2rad(2), deg2rad(-14), deg2rad(2), deg2rad(68), deg2rad(27.5), deg2rad(68)];
trafo_vetor = ["DiY", "DiY", "DiY", "YY", "YY", "YY"];

for kk=1:1
    %calcula as coisas
    trafoo = trafo_vetor(kk);
    l.pr.d = d_vetor(kk);
    l.pr.phi = phi_vetor(kk);

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
    [efc,ptot,Pm,cSw_p,cSw_s,sSw_p,sSw_s,Pil_cu,PiL_cu,PiLd_cu,Psc,P_core_L,P_core_tr] = f_loss(trafoo, out, l);
    Pm
    
    [out] = f_equations_plot(l,trafoo);
    C = num2cell(out);
    vPP(1:3,:) = cell2mat(C(1:3,:));
    vSS(1:3,:) = cell2mat(C(4:6,:));
    vLLdab(1:3,:) = cell2mat(C(7:9,:));
    hbb(1:3,:) = cell2mat(C(10:12,:));
    HBB(1:3,:) = cell2mat(C(13:15,:));
    idab_plot(1:3,:) = cell2mat(C(16:18,:));
    x0s(1:6,:) = cell2mat(C(19:24,:));
    ts(1,:) = cell2mat(C(25,:));
    sec_switch_plot(1,:) = cell2mat(C(26,:));
    sf_plot(1:6,:) = cell2mat(C(27:32,:));
    
    fase = 1;
    
    %% primario
    figure
    cmap = f_create_cmap(2, color2, color1);
    colormap(cmap)
    jetcustom = cmap;
    hold on
    stairs(ts*1e6,vPP(fase,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
    ylabel('$V_{L2_{A}}\,$[V]')
    yyaxis right
    plot(ts*1e6,x0s(fase,:),'Color',jetcustom(2,:), 'LineWidth',1.5)
    hold off
    grid on
    grid minor
    set(gca, 'FontSize', 20)
    xlabel('$t[\mu$s]')
    ylabel('$i_{trA}\,$[A]')
    % ylim([-3 3])
    yyaxis left
    % ylim([-300 300])
    ax = gca;
    ax.YAxis(1).Color = jetcustom(1,:);
    ax.YAxis(2).Color = jetcustom(2,:);
    title_name = append(trafoo,', primario, d=',string(l.pr.d),', phi = ',string(rad2deg(l.pr.phi)));
    title(title_name)
    file_name = append('figure\finalCap2\comparacoes\primario_exp_teorico_',trafoo,'_d_',string(l.pr.d),'_phi_',string(rad2deg(l.pr.phi)),'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
    
    %% secundario
    figure
    cmap = f_create_cmap(2, color2, color1);
    colormap(cmap)
    jetcustom = cmap;
    hold on
    stairs(ts*1e6,vSS(fase,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
    ylabel('$V_{L2_{A}}\,$[V]')
    yyaxis right
    plot(ts*1e6,x0s(3+fase,:),'Color',jetcustom(2,:), 'LineWidth',1.5)
    hold off
    grid on
    grid minor
    set(gca, 'FontSize', 20)
    xlabel('$t[\mu$s]')
    ylabel('$i_{trA}\,$[A]')
    % ylim([-3 3])
    yyaxis left
    % ylim([-300 300])
    ax = gca;
    ax.YAxis(1).Color = jetcustom(1,:);
    ax.YAxis(2).Color = jetcustom(2,:);
    title_name = append(trafoo,', secundario, d=',string(l.pr.d),', phi = ',string(rad2deg(l.pr.phi)));
    title(title_name)
    file_name = append('figure\finalCap2\comparacoes\secundario_exp_teorico_',trafoo,'_d_',string(l.pr.d),'_phi_',string(rad2deg(l.pr.phi)),'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
    %% indutor
    figure
    cmap = f_create_cmap(3, color2, color1);
    colormap(cmap)
    jetcustom = cmap;
    hold on
%     stairs(ts*1e6,vLLdab(fase,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
%     ylabel('$V_{L2_{A}}\,$[V]')
%     yyaxis right
    plot(ts*1e6,idab_plot(1,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
    plot(ts*1e6,idab_plot(2,:),'Color',jetcustom(2,:), 'LineWidth',1.5)
    plot(ts*1e6,idab_plot(3,:),'Color',jetcustom(3,:), 'LineWidth',1.5)
    hold off
    grid on
    grid minor
    set(gca, 'FontSize', 20)
    xlabel('$t[\mu$s]')
    ylabel('$i_{Ldab}\,$[A]')
    % ylim([-3 3])
    % ylim([-300 300])
    title_name = append(trafoo,', indutor, d=',string(l.pr.d),', phi = ',string(rad2deg(l.pr.phi)));
    title(title_name)
    file_name = append('figure\finalCap2\comparacoes\indutor_exp_teorico_',trafoo,'_d_',string(l.pr.d),'_phi_',string(rad2deg(l.pr.phi)),'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
    %% corrente das duas fases
    figure
    cmap = f_create_cmap(2, color2, color1);
    colormap(cmap)
    jetcustom = cmap;
    hold on
    plot(ts*1e6,x0s(fase,:),'Color',jetcustom(1,:), 'LineWidth',1.5)
    plot(ts*1e6,-x0s(3+fase,:),'Color',jetcustom(2,:), 'LineWidth',1.5)
    hold off
    grid on
    grid minor
    set(gca, 'FontSize', 20)
    xlabel('$t[\mu$s]')
    ylabel('$i_{trA}\,$[A]')
    ax = gca;
    title_name = append(trafoo,', prim e -sec, d=',string(l.pr.d),', phi = ',string(rad2deg(l.pr.phi)));
    title(title_name)
    file_name = append('figure\finalCap2\comparacoes\pri_menos_sec_exp_teorico_',trafoo,'_d_',string(l.pr.d),'_phi_',string(rad2deg(l.pr.phi)),'.pdf');
    exportgraphics(gca,file_name,'ContentType','vector');
end








