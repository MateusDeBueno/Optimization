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


load('YY.mat');
load('YD.mat');
load('DfD.mat');
load('DiD.mat');
load('DiY.mat');
load('DfY.mat');
% Combine the structs into a single struct
l.eq = struct('YY', YY, 'YD', YD, 'DfD', DfD, ...
              'DiD', DiD, 'DiY', DiY, 'DfY', DfY);

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

l.pr.phi = deg2rad(15);


out = f_equationsDfD_dt(l);
C = num2cell(out);
[hbrm,HBrm,Ip,Is,iiRMS,iiME,ioRMS,ioME,Pm,idrm,ilrm,iLrm,iSwPrm,iSwSrm,P_core_tr,Bpk_tr,P_core_L,Bpk_L] = C{:}
