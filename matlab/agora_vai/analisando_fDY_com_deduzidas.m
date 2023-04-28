clear; 
close all; 
clc;

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
%%

trafo = 'oDY';

color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue

%% analisar teorico com as equacoes deduzidas




l.L.a = 1.285676564102240;
l.L.b = 2.349349709620920;
l.L.kc = 11.187048860336922;
l.L.int_ki = integral(@(theta) abs(cos(theta)).^l.L.a,0,2*pi);
l.L.ki = l.L.kc/(2^(l.L.b-l.L.a)*(2*pi)^(l.L.a-1)*l.L.int_ki);
l.L.Ac = 421.3e-6;
l.L.Ve = 52.1e-6;
l.L.N = 14;
l.L.strand = 180;
l.L.awg = 38;
l.L.dl = 0; % [m]
l.L.MLT = 140*1e-3; % [m]
l.L.Nl = 1;

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
l.tr.dl = 0; % [m]
l.tr.MLT = 230*1e-3; % [m]

l.sC.Rac100 = (4.3e-3)/2;
l.C.Rb_ac10k = 4.2e-3 + 2e-3;
l.pr.dt = 0e-9;
l.pr.Vi = 400;
l.pr.d = 0.75;
l.pr.fs = 100e3;
l.pr.Ldab = 60e-6;
l.pr.Ld1 = 2e-6;
l.pr.n = 1;
l.pr.Ld2 = l.pr.Ld1*l.pr.n*l.pr.n;
l.pr.Lm = 500e-6;
l.pr.M = l.pr.Lm*l.pr.n;
l.pr.L1 = l.pr.Ld1 + l.pr.Lm;
l.pr.L2 = l.pr.Ld2 + l.pr.n*l.pr.n*l.pr.Lm;
l.pr.k = l.pr.M/sqrt(l.pr.L1.*l.pr.L2);

l.pr.phi = 0;

load('YY.mat');
load('YD.mat');
load('DfD.mat');
load('DiD.mat');
load('DiY.mat');
load('DfY.mat');
% Combine the structs into a single struct
l.eq = struct('YY', YY, 'YD', YD, 'DfD', DfD, ...
              'DiD', DiD, 'DiY', DiY, 'DfY', DfY);

pr = 360*2; %precision
vec_ph = -180:360/pr:180;

for kk=1:length(vec_ph)
    l.pr.phi = deg2rad(vec_ph(kk));
    out = f_equationsDfY_dt(l);
    C = num2cell(out);
    [hbrm,HBrm,Ip,Is,iiRMS,iiME(kk),ioRMS(kk),ioME(kk),Pm(kk),idrm(kk),ilrm(kk),iLrm(kk),iSwPrm,iSwSrm,P_core_tr,Bpk_tr,P_core_L,Bpk_L] = C{:};
end

%%


figure

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

hold on

x_points = [-29.9, -29.9, 60, 60];  
y_points = [-1500, 1500, 1500, -1500];
color = [0, 0, 1];
a = fill(x_points, y_points, color, 'LineWidth',.06);
a.FaceAlpha = 0.08;
a.EdgeAlpha = 0.4;
text(60,-2000/2,'$\leftarrow$ Region 1','Interpreter', 'Latex','FontSize', 14) 

x_points = [-120, -120, -30.1, -30.1];  
y_points = [-1500, 1500, 1500, -1500];
color = [0, 0, 1];
a = fill(x_points, y_points, color, 'LineWidth',.06);
a.FaceAlpha = 0.08;
a.EdgeAlpha = 0.4;
text(-124,-500/2,'Region','Interpreter', 'Latex','FontSize', 14,'HorizontalAlignment', 'right') 
text(-124,-1000/2,'2 $\rightarrow$','Interpreter', 'Latex','FontSize', 14,'HorizontalAlignment', 'right') 


plot(vec_ph,Pm, 'Color',jetcustom(1,:), 'LineWidth',1.5)

ylabel('$P_o$[W]')
ylim([-1200 1200])
yticks(-1200:600:1200)

yyaxis right
ylabel('$i_{Ldab-rms}\,[A]$')
plot(vec_ph,idrm, 'Color',jetcustom(2,:), 'LineWidth',1.5)
ylim([2 8])
yticks(2:2:8)
hold off
xlabel('$\phi\,[^{\circ}]$')
set(gca, 'FontSize', 20)
grid on
grid minor
ax = gca;
ax.YAxis(1).Color = jetcustom(1,:);
ax.YAxis(2).Color = jetcustom(2,:);
xlim([-180 180])
xticks(-180:60:180)

file_name = append('figure\finalCap2\Po_bid_',trafo,'.pdf');
exportgraphics(gca,file_name,'ContentType','vector');
