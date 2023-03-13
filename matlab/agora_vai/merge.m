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

color1 = [0.045 0.245 0.745]; % blue
color2 = [0.635 0.635 0.635]; % gray

syms Ld1 Ld2 n Lm Po t L1 L2 Ldab M fs Vi d dt real positive
syms phi real

Vo = Vi*d;
Ts = 1/fs;

phi_num = [deg2rad(-10),deg2rad(50)];  %[MUDAR]
pr = 200; %precision

%cria intervalo de angulo
YD.intervalo.f1 = [deg2rad(-30) 0];
YD.intervalo.f2 = [0 deg2rad(60)];
[YD.x0s.f1, YD.ts.f1, YD.idab.f1, YD.hb.f1, YD.HB.f1, YD.Ip.f1, YD.Is.f1, YD.iME.f1, YD.idrm.f1, YD.ilrm.f1, YD.iLrm.f1, YD.iSwPrm.f1, YD.iSwSrm.f1] = simplify_YD(min((YD.intervalo.f1+0.1)));
[YD.x0s.f2, YD.ts.f2, YD.idab.f2, YD.hb.f2, YD.HB.f2, YD.Ip.f2, YD.Is.f2, YD.iME.f2, YD.idrm.f2, YD.ilrm.f2, YD.iLrm.f2, YD.iSwPrm.f2, YD.iSwSrm.f2] = simplify_YD(min((YD.intervalo.f2+0.1)));
YD.vec_ph.f2 = min(YD.intervalo.f2):(max(YD.intervalo.f2)-min(YD.intervalo.f2))/pr:max(YD.intervalo.f2);
YD.vec_ph.f1 = min(YD.intervalo.f1):(max(YD.intervalo.f1)-min(YD.intervalo.f1))/pr:max(YD.intervalo.f1);

%equacao_limite
YD.Ip_eq.f1 = YD.Ip.f1 == 0;
YD.Ip_eq.f2 = YD.Ip.f2 == 0;
YD.Is_eq.f1 = YD.Is.f1 == 0;
YD.Is_eq.f2 = YD.Is.f2 == 0;
YD.pot_eq.f1 = YD.iME.f1*Vo == Po;
YD.pot_eq.f2 = YD.iME.f2*Vo == Po;

%cria funcao
YD.fpot_eq.f1 = matlabFunction(solve(YD.pot_eq.f1,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIp_eq.f1 = matlabFunction(solve(YD.Ip_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIs_eq.f1 = matlabFunction(solve(YD.Is_eq.f1,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fpot_eq.f2 = matlabFunction(solve(YD.pot_eq.f2,Po), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIp_eq.f2 = matlabFunction(solve(YD.Ip_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
YD.fIs_eq.f2 = matlabFunction(solve(YD.Is_eq.f2,d), 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});


%% plot Po por phi

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
hold on
%legend that is unrelated to the plotted data
L1 = plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
L2 = plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
L3 = plot(nan, nan,'Color',jetcustom(3,:),'LineWidth',1.5);

lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),lim_pot.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_pot.f2,'Color',jetcustom(1,:),'LineWidth',1.5)

lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),lim_pot.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_pot.f2,'Color',jetcustom(2,:),'LineWidth',1.5)


lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),lim_pot.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_pot.f2,'Color',jetcustom(3,:),'LineWidth',1.5)

hold off
ylim([0 1500])
xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlabel('$\phi\,[^{\circ}]$')
ylabel('$P_o\,$[W]')
set(gca, 'FontSize', 20)
