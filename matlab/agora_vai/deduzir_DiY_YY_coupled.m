clear; close all; clc;

YY.trafo = 'YY';
DiY.trafo = 'DiY';


% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

addpath('utils')
addpath('utils_transf')

% color1 = [0.045 0.245 0.745]; % blue
% color2 = [0.635 0.635 0.635]; % gray
color1 = [249,152,32]/255; % orange
color2 = [32,129,249]/255; % blue


syms Ld1 Ld2 n Lm Po t L1 L2 Ldab M fs Vi d dt real positive
syms phi real

Vo = Vi*d;
Ts = 1/fs;

YY.phi_num = [deg2rad(10),deg2rad(80)];  %[MUDAR]
[YY.x0s.f1, YY.ts.f1, YY.idab.f1, YY.hb.f1, YY.HB.f1, YY.Ip.f1, YY.Is.f1, YY.iME.f1, YY.idrm.f1, YY.ilrm.f1, YY.iLrm.f1, YY.iSwPrm.f1, YY.iSwSrm.f1] = simplify_YY(YY.phi_num(1));
[YY.x0s.f2, YY.ts.f2, YY.idab.f2, YY.hb.f2, YY.HB.f2, YY.Ip.f2, YY.Is.f2, YY.iME.f2, YY.idrm.f2, YY.ilrm.f2, YY.iLrm.f2, YY.iSwPrm.f2, YY.iSwSrm.f2] = simplify_YY(YY.phi_num(2));

DiY.phi_num = [deg2rad(-10),deg2rad(50)];  %[MUDAR]
[DiY.x0s.f1, DiY.ts.f1, DiY.idab.f1, DiY.hb.f1, DiY.HB.f1, DiY.Ip.f1, DiY.Is.f1, DiY.iME.f1, DiY.idrm.f1, DiY.ilrm.f1, DiY.iLrm.f1, DiY.iSwPrm.f1, DiY.iSwSrm.f1] = simplify_DiY(DiY.phi_num(1));
[DiY.x0s.f2, DiY.ts.f2, DiY.idab.f2, DiY.hb.f2, DiY.HB.f2, DiY.Ip.f2, DiY.Is.f2, DiY.iME.f2, DiY.idrm.f2, DiY.ilrm.f2, DiY.iLrm.f2, DiY.iSwPrm.f2, DiY.iSwSrm.f2] = simplify_DiY(DiY.phi_num(2));


Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61.5e-6;
Ld1_num = 1.4e-6;
n_num = 5/9;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = [350e-6, .7e-3 10e-3];
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);


%%
%cria intervalo de angulo
YY.intervalo.f1 = [0 deg2rad(60)];
YY.intervalo.f2 = [deg2rad(60) deg2rad(89.999999)];

%% passo para plot
pr = 200; %precision
YY.vec_ph.f2 = min(YY.intervalo.f2):(max(YY.intervalo.f2)-min(YY.intervalo.f2))/pr:max(YY.intervalo.f2);
YY.vec_ph.f1 = min(YY.intervalo.f1):(max(YY.intervalo.f1)-min(YY.intervalo.f1))/pr:max(YY.intervalo.f1);

YY.fx0s.f1 = matlabFunction(YY.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fx0s.f2 = matlabFunction(YY.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fts.f1 = matlabFunction(YY.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fts.f2 = matlabFunction(YY.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%equacao_limite
YY.Ip_eq.f1 = YY.Ip.f1 == 0;
YY.Is_eq.f1 = YY.Is.f1 == 0;
YY.pot_eq.f1 = YY.iME.f1*Vo == Po;

%%

YY.fpot_eq.f1 = matlabFunction(solve(YY.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIp_eq.f1 = matlabFunction(solve(YY.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIs_eq.f1 = matlabFunction(solve(YY.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

YY.Ip_eq.f2 = YY.Ip.f2 == 0;
YY.Is_eq.f2 = YY.Is.f2 == 0;
YY.pot_eq.f2 = YY.iME.f2*Vo == Po;

YY.fpot_eq.f2 = matlabFunction(solve(YY.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIp_eq.f2 = matlabFunction(solve(YY.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIs_eq.f2 = matlabFunction(solve(YY.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%correntes eficazes
YY.fidrm.f1 = matlabFunction(YY.idrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.filrm.f1 = matlabFunction(YY.ilrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fiLrm.f1 = matlabFunction(YY.iLrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fiSwPrm.f1 = matlabFunction(YY.iSwPrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fiSwSrm.f1 = matlabFunction(YY.iSwSrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

YY.fidrm.f2 = matlabFunction(YY.idrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.filrm.f2 = matlabFunction(YY.ilrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fiLrm.f2 = matlabFunction(YY.iLrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fiSwPrm.f2 = matlabFunction(YY.iSwPrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fiSwSrm.f2 = matlabFunction(YY.iSwSrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

YY.fIp.f1 = matlabFunction(YY.Ip.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIs.f1 = matlabFunction(YY.Is.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIp.f2 = matlabFunction(YY.Ip.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YY.fIs.f2 = matlabFunction(YY.Is.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

mag=1;

YY.fIp.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.phi_num(1))
YY.fIs.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.phi_num(1))
YY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.phi_num(1))

YY.fIp.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.phi_num(2))
YY.fIs.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.phi_num(2))
YY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.phi_num(2))

%%
%cria intervalo de angulo
DiY.intervalo.f1 = [deg2rad(-30) 0];
DiY.intervalo.f2 = [0 deg2rad(59.9999999)];

%% passo para plot
pr = 200; %precision
DiY.vec_ph.f2 = min(DiY.intervalo.f2):(max(DiY.intervalo.f2)-min(DiY.intervalo.f2))/pr:max(DiY.intervalo.f2);
DiY.vec_ph.f1 = min(DiY.intervalo.f1):(max(DiY.intervalo.f1)-min(DiY.intervalo.f1))/pr:max(DiY.intervalo.f1);

DiY.fx0s.f1 = matlabFunction(DiY.x0s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fx0s.f2 = matlabFunction(DiY.x0s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fts.f1 = matlabFunction(DiY.ts.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fts.f2 = matlabFunction(DiY.ts.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

%equacao_limite
DiY.Ip_eq.f1 = DiY.Ip.f1 == 0;
DiY.Is_eq.f1 = DiY.Is.f1 == 0;
DiY.pot_eq.f1 = DiY.iME.f1*Vo == Po;

%%

DiY.fpot_eq.f1 = matlabFunction(solve(DiY.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fIp_eq.f1 = matlabFunction(solve(DiY.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fIs_eq.f1 = matlabFunction(solve(DiY.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

DiY.Ip_eq.f2 = DiY.Ip.f2 == 0;
DiY.Is_eq.f2 = DiY.Is.f2 == 0;
DiY.pot_eq.f2 = DiY.iME.f2*Vo == Po;

DiY.fpot_eq.f2 = matlabFunction(solve(DiY.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fIp_eq.f2 = matlabFunction(solve(DiY.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fIs_eq.f2 = matlabFunction(solve(DiY.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%correntes eficazes
DiY.fidrm.f1 = matlabFunction(DiY.idrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.filrm.f1 = matlabFunction(DiY.ilrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fiLrm.f1 = matlabFunction(DiY.iLrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fiSwPrm.f1 = matlabFunction(DiY.iSwPrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fiSwSrm.f1 = matlabFunction(DiY.iSwSrm.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

DiY.fidrm.f2 = matlabFunction(DiY.idrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.filrm.f2 = matlabFunction(DiY.ilrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fiLrm.f2 = matlabFunction(DiY.iLrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fiSwPrm.f2 = matlabFunction(DiY.iSwPrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fiSwSrm.f2 = matlabFunction(DiY.iSwSrm.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

DiY.fIp.f1 = matlabFunction(DiY.Ip.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fIs.f1 = matlabFunction(DiY.Is.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fIp.f2 = matlabFunction(DiY.Ip.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
DiY.fIs.f2 = matlabFunction(DiY.Is.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});

mag=1;

DiY.fIp.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.phi_num(1))
DiY.fIs.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.phi_num(1))
DiY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.phi_num(1))

DiY.fIp.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.phi_num(2))
DiY.fIs.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.phi_num(2))
DiY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.phi_num(2))




%% plot d por Po
% 
% cmap = f_create_cmap(2, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% 
% figure
% hold on
% 
% L11 = plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
% L22 = plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
% 
% YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
% YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
% YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
% 
% fill([YY.lim_pot.f1 fliplr(YY.lim_pot.f1)], [YY.lim_p.f1 fliplr(YY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
% YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
% YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
% YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
% 
% fill([YY.lim_pot.f2 fliplr(YY.lim_pot.f2)], [YY.lim_p.f2 fliplr(YY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
% plot(YY.lim_pot.f1,YY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(YY.lim_pot.f2,YY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(YY.lim_pot.f1,YY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(YY.lim_pot.f2,YY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% 
% DiY.lim_p.f1 = ones(1,length(DiY.vec_ph.f1)).*DiY.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f1);
% DiY.lim_s.f1 = ones(1,length(DiY.vec_ph.f1)).*DiY.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f1);
% DiY.lim_pot.f1 = ones(1,length(DiY.vec_ph.f1)).*DiY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f1);
% 
% fill([DiY.lim_pot.f1 fliplr(DiY.lim_pot.f1)], [DiY.lim_p.f1 fliplr(DiY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
% DiY.lim_p.f2 = ones(1,length(DiY.vec_ph.f2)).*DiY.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f2);
% DiY.lim_s.f2 = ones(1,length(DiY.vec_ph.f2)).*DiY.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f2);
% DiY.lim_pot.f2 = ones(1,length(DiY.vec_ph.f2)).*DiY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f2);
% 
% fill([DiY.lim_pot.f2 fliplr(DiY.lim_pot.f2)], [DiY.lim_p.f2 fliplr(DiY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
% plot(DiY.lim_pot.f1,DiY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
% plot(DiY.lim_pot.f2,DiY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
% plot(DiY.lim_pot.f1,DiY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
% plot(DiY.lim_pot.f2,DiY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
% 
% hold off
% ylim([0 1.4])
% xlim([0 max(YY.lim_pot.f2)])
% legend({'YY','DiY'},'Location','best','FontSize', 14)
% 
% grid on
% grid minor
% xlabel('$P_o$[W]')
% ylabel('$d$')
% set(gca, 'FontSize', 20)


%% plot Vi por Po


% mag = 3;
% for mag=1:1:3
for mag=1:1:3
    figure
    hold on
    
    cmap = f_create_cmap(2, color2, color1);
    colormap(cmap)
    jetcustom = cmap;
    
    L11 = plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
    L22 = plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
    
    YY.lim_p.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
    YY.lim_s.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
    YY.lim_pot.f1 = ones(1,length(YY.vec_ph.f1)).*YY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f1);
    
    fill([YY.lim_pot.f1 fliplr(YY.lim_pot.f1)], [Vi_num*YY.lim_p.f1 fliplr(Vi_num*YY.lim_s.f1)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
    
    YY.lim_p.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
    YY.lim_s.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
    YY.lim_pot.f2 = ones(1,length(YY.vec_ph.f2)).*YY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,YY.vec_ph.f2);
    
    fill([YY.lim_pot.f2 fliplr(YY.lim_pot.f2)], [Vi_num*YY.lim_p.f2 fliplr(Vi_num*YY.lim_s.f2)],jetcustom(1,:), 'FaceAlpha', .3, 'EdgeColor', 'none');
    
    plot(YY.lim_pot.f1,Vi_num*YY.lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
    plot(YY.lim_pot.f2,Vi_num*YY.lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
    plot(YY.lim_pot.f1,Vi_num*YY.lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
    plot(YY.lim_pot.f2,Vi_num*YY.lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
    
    DiY.lim_p.f1 = ones(1,length(DiY.vec_ph.f1)).*DiY.fIp_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f1);
    DiY.lim_s.f1 = ones(1,length(DiY.vec_ph.f1)).*DiY.fIs_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f1);
    DiY.lim_pot.f1 = ones(1,length(DiY.vec_ph.f1)).*DiY.fpot_eq.f1(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f1);
    
    fill([DiY.lim_pot.f1 fliplr(DiY.lim_pot.f1)], [Vi_num*DiY.lim_p.f1 fliplr(Vi_num*DiY.lim_s.f1)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
    
    DiY.lim_p.f2 = ones(1,length(DiY.vec_ph.f2)).*DiY.fIp_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f2);
    DiY.lim_s.f2 = ones(1,length(DiY.vec_ph.f2)).*DiY.fIs_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f2);
    DiY.lim_pot.f2 = ones(1,length(DiY.vec_ph.f2)).*DiY.fpot_eq.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,DiY.vec_ph.f2);
    
    fill([DiY.lim_pot.f2 fliplr(DiY.lim_pot.f2)], [Vi_num*DiY.lim_p.f2 fliplr(Vi_num*DiY.lim_s.f2)],jetcustom(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
    
    plot(DiY.lim_pot.f1,Vi_num*DiY.lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
    plot(DiY.lim_pot.f2,Vi_num*DiY.lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
    plot(DiY.lim_pot.f1,Vi_num*DiY.lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
    plot(DiY.lim_pot.f2,Vi_num*DiY.lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
    
    hold off
    ylim([0 1.4*Vi_num])
    xlim([0 max(YY.lim_pot.f2)])
    legend({'YY','DiY'},'Location','best','FontSize', 14)
    
    grid on
    grid minor
    xlabel('$P_o$[W]')
    ylabel('$V_o[V]$')
    set(gca, 'FontSize', 20)
    title('ZVS range')
end