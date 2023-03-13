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

phi_num = [deg2rad(10),deg2rad(80)];  %[MUDAR]
[x0s.f1, ts.f1, idab.f1, hb.f1, HB.f1, Ip.f1, Is.f1, iME.f1, idrm.f1, ilrm.f1, iLrm.f1, iSwPrm.f1, iSwSrm.f1] = simplify_YY(phi_num(1));
[x0s.f2, ts.f2, idab.f2, hb.f2, HB.f2, Ip.f2, Is.f2, iME.f2, idrm.f2, ilrm.f2, iLrm.f2, iSwPrm.f2, iSwSrm.f2] = simplify_YY(phi_num(2));

trafo = 'YY';

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
n_num = 1;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = [700e-6, 1.4e-3 10e-3];
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

%cria intervalo de angulo
intervalo.f1 = [0 deg2rad(60)];
intervalo.f2 = [deg2rad(60) deg2rad(89.999999)];

%% passo para plot
pr = 200; %precision
vec_ph.f2 = min(intervalo.f2):(max(intervalo.f2)-min(intervalo.f2))/pr:max(intervalo.f2);
vec_ph.f1 = min(intervalo.f1):(max(intervalo.f1)-min(intervalo.f1))/pr:max(intervalo.f1);

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

%%
for mag=1:3
    equat.f2 = fx0s.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,max(intervalo.f1));
    equat_ts.f2 = fts.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,max(intervalo.f1));
    
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
    legend({'$i_l$','$-i_L$'},'Location','best','FontSize', 14)
    set(gca, 'FontSize', 20)
    xlabel('$t[\mu$s]')
    ylabel('$i\,$[A]')
    f_save_figure(append('figure\',string(trafo),'_',string(mag),'_ilIL.pdf'))
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

%%

% Create the contour plot
figure
hold on
[C,h] = contour(x1*1e3, y1, z1,[0.94 0.96 0.98 1 1.005 1.01 1.04 1.1], 'ShowText', 'on','Color',color1);
h.LineWidth = 1.5;
clabel(C,h,'FontSize',15,'Color',[0 0 0],'LineWidth',2)
[C,h] = contour(x2*1e3, y2, z2,[0.94 0.96 0.98 1 1.005 1.01 1.04 1.1], 'ShowText', 'on','Color',color1);
h.LineWidth = 1.5;
clabel(C,h,'FontSize',15,'Color',[0 0 0],'LineWidth',2)
hold off
grid on
grid minor
xlabel('$L_m\,$[mH]')
ylabel('$\phi\,[^{\circ}]$')
set(gca, 'FontSize', 20)
f_save_figure(append('figure\',string(trafo),'_LM_ratio.pdf'))
%% plot d por phi
figure
cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;
hold on
%legend that is unrelated to the plotted data
L1 = plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
L2 = plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
L3 = plot(nan, nan,'Color',jetcustom(3,:),'LineWidth',1.5);

text(13,1.7,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(13,0.35,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(60,1,'SS in both','Interpreter', 'Latex','FontSize', 14) 

lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f1),lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)

fill([rad2deg(vec_ph.f1) fliplr(rad2deg(vec_ph.f1))], [lim_p.f1 fliplr(lim_s.f1)],jetcustom(1,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([rad2deg(vec_ph.f2) fliplr(rad2deg(vec_ph.f2))], [lim_p.f2 fliplr(lim_s.f2)],jetcustom(1,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f1),lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)

fill([rad2deg(vec_ph.f1) fliplr(rad2deg(vec_ph.f1))], [lim_p.f1 fliplr(lim_s.f1)],jetcustom(2,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([rad2deg(vec_ph.f2) fliplr(rad2deg(vec_ph.f2))], [lim_p.f2 fliplr(lim_s.f2)],jetcustom(2,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),lim_p.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_p.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f1),lim_s.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),lim_s.f2,'Color',jetcustom(3,:),'LineWidth',1.5)

fill([rad2deg(vec_ph.f1) fliplr(rad2deg(vec_ph.f1))], [lim_p.f1 fliplr(lim_s.f1)],jetcustom(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([rad2deg(vec_ph.f2) fliplr(rad2deg(vec_ph.f2))], [lim_p.f2 fliplr(lim_s.f2)],jetcustom(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

hold off
ylim([0 3])
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlabel('$\phi\,[^{\circ}]$')
ylabel('$d$')
set(gca, 'FontSize', 20)
xlim([min(intervalo.f1)*180/pi max(intervalo.f2)*180/pi])
f_save_figure(append('figure\',string(trafo),'_d_phi.pdf'))
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
ylim([0 2500])
xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlabel('$\phi\,[^{\circ}]$')
ylabel('$P_o\,$[W]')
set(gca, 'FontSize', 20)
f_save_figure(append('figure\',string(trafo),'_gain_static.pdf'))

%% plot d por Po

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
hold on
%legend that is unrelated to the plotted data
L1 = plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
L2 = plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
L3 = plot(nan, nan,'Color',jetcustom(3,:),'LineWidth',1.5);

text(500,1.7,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(500,0.35,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(1500,1,'SS in both','Interpreter', 'Latex','FontSize', 14) 

lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);

fill([lim_pot.f1 fliplr(lim_pot.f1)], [lim_p.f1 fliplr(lim_s.f1)],jetcustom(1,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);

fill([lim_pot.f2 fliplr(lim_pot.f2)], [lim_p.f2 fliplr(lim_s.f2)],jetcustom(1,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

plot(lim_pot.f1,lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(lim_pot.f2,lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(lim_pot.f1,lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(lim_pot.f2,lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)


lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);

fill([lim_pot.f1 fliplr(lim_pot.f1)], [lim_p.f1 fliplr(lim_s.f1)],jetcustom(2,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);

fill([lim_pot.f2 fliplr(lim_pot.f2)], [lim_p.f2 fliplr(lim_s.f2)],jetcustom(2,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

plot(lim_pot.f1,lim_p.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(lim_pot.f2,lim_p.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(lim_pot.f1,lim_s.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(lim_pot.f2,lim_s.f2,'Color',jetcustom(2,:),'LineWidth',1.5)

lim_p.f1 = ones(1,length(vec_ph.f1)).*fIp_eq.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
lim_s.f1 = ones(1,length(vec_ph.f1)).*fIs_eq.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
lim_pot.f1 = ones(1,length(vec_ph.f1)).*fpot_eq.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);

fill([lim_pot.f1 fliplr(lim_pot.f1)], [lim_p.f1 fliplr(lim_s.f1)],jetcustom(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

lim_p.f2 = ones(1,length(vec_ph.f2)).*fIp_eq.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
lim_s.f2 = ones(1,length(vec_ph.f2)).*fIs_eq.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
lim_pot.f2 = ones(1,length(vec_ph.f2)).*fpot_eq.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);

fill([lim_pot.f2 fliplr(lim_pot.f2)], [lim_p.f2 fliplr(lim_s.f2)],jetcustom(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

plot(lim_pot.f1,lim_p.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(lim_pot.f2,lim_p.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(lim_pot.f1,lim_s.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(lim_pot.f2,lim_s.f2,'Color',jetcustom(3,:),'LineWidth',1.5)

hold off
ylim([0 3])
xlim([0 max(lim_pot.f2)])
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlabel('$P_o$[W]')
ylabel('$d$')
set(gca, 'FontSize', 20)
f_save_figure(append('figure\',string(trafo),'_d_Po.pdf'))
