clear; close all; clc;

name = 'DfY';

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

syms fs Vi d Po dt t real positive
syms a1 b1 z1 a2 b2 z2 real
syms xa1 xb1 xz1 xa2 xb2 xz2 real
syms phi real
syms Ld1 Ld2 n Lm real positive
syms L1 L2 Ldab M real positive
syms fs Vi d dt real positive
syms phi real

Vo = Vi*d;
Ts = 1/fs;

Vi_num = 400;
d_num = 1;
fs_num = 100e3;
Ldab_num = 61e-6;
Ld1_num = 2e-6;
% Ld1_num = 0;
n_num = 1;
Ld2_num = Ld1_num*n_num*n_num;
Lm_num = [700e-6, 1.4e-3 10e-3];
M_num = Lm_num*n_num;
L1_num = Ld1_num + Lm_num;
L2_num = Ld2_num + n_num*n_num*Lm_num;
k_num = M_num/sqrt(L1_num.*L2_num);

phi_num = [deg2rad(10),deg2rad(80)];  %[MUDAR]
%cria intervalo de angulo
intervalo.f1 = [0 deg2rad(60)];
intervalo.f2 = [deg2rad(60) deg2rad(89.999999)];



[x0s.f1, ts.f1, idab.f1, hb.f1, HB.f1, Ip.f1, Is.f1, iME.f1, idrm.f1, ilrm.f1, iLrm.f1, iSwPrm.f1, iSwSrm.f1,B.f1] = simplify_YY(phi_num(1));
[x0s.f2, ts.f2, idab.f2, hb.f2, HB.f2, Ip.f2, Is.f2, iME.f2, idrm.f2, ilrm.f2, iLrm.f2, iSwPrm.f2, iSwSrm.f2,B.f2] = simplify_YY(phi_num(2));



%%



%passo para plot
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

mag=1;

equat.f2 = fx0s.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(2));
equat_ts.f2 = fts.f2(L1_num(mag),L2_num(mag),Ldab_num,M_num(mag),Vi_num,d_num,fs_num,phi_num(2));

cmap = f_create_cmap(2, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
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
f_save_figure(append('figure\YY_i',string(mag),'_ilIL.pdf'))


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
f_save_figure('figure\YY_LM_ratio.pdf')

%% plot d por phi

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
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
xlim([min(intervalo.f1)*180/pi max(intervalo.f2)]*180/pi)
f_save_figure('figure\YY_d_phi.pdf')




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
f_save_figure('figure\YY_gain_static.pdf')



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

text(300,2.5,'HS in primary','Interpreter', 'Latex','FontSize', 14) 
text(800,0.35,'HS in secondary','Interpreter', 'Latex','FontSize', 14) 
text(1000,1.55,'SS in both','Interpreter', 'Latex','FontSize', 14) 

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
f_save_figure('figure\YY_d_Po.pdf')




%% plot iLdab por phi
% fidrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,deg2rad(60))
% filrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,deg2rad(60))
% fiLrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,deg2rad(60))
% fiSwPrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,deg2rad(60))
% fiSwSrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,deg2rad(60))

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
hold on
%legend that is unrelated to the plotted data
L1 = plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
L2 = plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
L3 = plot(nan, nan,'Color',jetcustom(3,:),'LineWidth',1.5);

equat.f1 = ones(1,length(vec_ph.f1)).*fidrm.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*fidrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(1,:),'LineWidth',1.5)

equat.f1 = ones(1,length(vec_ph.f1)).*fidrm.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*fidrm.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(2,:),'LineWidth',1.5)

equat.f1 = ones(1,length(vec_ph.f1)).*fidrm.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*fidrm.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(3,:),'LineWidth',1.5)

hold off
xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
grid on
grid minor
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
xlabel('$\phi\,[^{\circ}]$')
ylabel('$i_{rms}\,$[A]')
set(gca, 'FontSize', 20)
f_save_figure('figure\YY_ILdab_phi.pdf')

%% corrente do primario e do secundario

cmap = f_create_cmap(3, color2, color1);
colormap(cmap)
jetcustom = cmap;

figure
hold on
%legend that is unrelated to the plotted data
% plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
% plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
% 
% plot(nan, nan,'Color',jetcustom(3,:),'LineWidth',1.5);
% plot(nan, nan,'Color',[0,0,0],'LineWidth',1.5);
% plot(nan, nan,'-.','Color',[0,0,0],'LineWidth',1.5);

plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);

plot(nan, nan,'Color',jetcustom(3,:),'LineWidth',1.5);

% plot(nan, nan,'-.','Color',jetcustom(1,:),'LineWidth',1.5);
% plot(nan, nan,'-.','Color',jetcustom(2,:),'LineWidth',1.5);
% 
% plot(nan, nan,'-.','Color',jetcustom(3,:),'LineWidth',1.5);
plot(nan, nan,'Color',[0,0,0],'LineWidth',1.5);
plot(nan, nan,'-.','Color',[0,0,0],'LineWidth',1.5);


%primario
equat.f1 = ones(1,length(vec_ph.f1)).*filrm.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*filrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(1,:),'LineWidth',1.5)

equat.f1 = ones(1,length(vec_ph.f1)).*filrm.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*filrm.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(2,:),'LineWidth',1.5)

equat.f1 = ones(1,length(vec_ph.f1)).*filrm.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*filrm.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(3,:),'LineWidth',1.5)

x1 = (vec_ph.f2(end)+vec_ph.f2(1))/2;
y1 = (equat.f2(end)+equat.f2(1))/2;
text(38,2.6,'$\leftarrow\,I_p$','Interpreter', 'Latex','FontSize', 14) 

%secundario
equat.f1 = ones(1,length(vec_ph.f1)).*fiLrm.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*fiLrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1, '-.','Color',jetcustom(1,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2, '-.','Color',jetcustom(1,:),'LineWidth',1.5)

equat.f1 = ones(1,length(vec_ph.f1)).*fiLrm.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*fiLrm.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1, '-.','Color',jetcustom(2,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2, '-.','Color',jetcustom(2,:),'LineWidth',1.5)

equat.f1 = ones(1,length(vec_ph.f1)).*fiLrm.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
equat.f2 = ones(1,length(vec_ph.f2)).*fiLrm.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
plot(rad2deg(vec_ph.f1),equat.f1, '-.','Color',jetcustom(3,:),'LineWidth',1.5)
plot(rad2deg(vec_ph.f2),equat.f2, '-.','Color',jetcustom(3,:),'LineWidth',1.5)

% x11 = (vec_ph.f2(end)+vec_ph.f2(1))/2;
% y11 = (equat.f2(end)+equat.f2(1))/2;
text(-4,4,'$I_s\,\rightarrow$','Interpreter', 'Latex','FontSize', 14) 

hold off
legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10\,$mH'},'Location','best','FontSize', 14,'NumColumns',1);
ylim([0.5 3.5])
xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
grid on
grid minor
xlabel('$\phi\,[^{\circ}]$')
ylabel('$i_{rms}\,$[A]')
set(gca, 'FontSize', 20)
f_save_figure('figure\YY_iLl_phi.pdf')


%%








%% plot da corrente do secundario refletida pro primario

% cmap = f_create_cmap(3, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% 
% figure
% hold on
% %legend that is unrelated to the plotted data
% plot(nan, nan,'Color',jetcustom(1,:),'LineWidth',1.5);
% plot(nan, nan,'Color',jetcustom(2,:),'LineWidth',1.5);
% 
% plot(nan, nan,'Color',jetcustom(3,:),'LineWidth',1.5);
% plot(nan, nan,'Color',[0,0,0],'LineWidth',1.5);
% plot(nan, nan,'-.','Color',[0,0,0],'LineWidth',1.5);
% 
% 
% 
% %primario
% equat.f1 = ones(1,length(vec_ph.f1)).*filrm.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
% equat.f2 = ones(1,length(vec_ph.f2)).*filrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
% plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% 
% equat.f1 = ones(1,length(vec_ph.f1)).*filrm.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
% equat.f2 = ones(1,length(vec_ph.f2)).*filrm.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
% plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
% plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
% 
% equat.f1 = ones(1,length(vec_ph.f1)).*filrm.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
% equat.f2 = ones(1,length(vec_ph.f2)).*filrm.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
% plot(rad2deg(vec_ph.f1),equat.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
% plot(rad2deg(vec_ph.f2),equat.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
% 
% x1 = (vec_ph.f2(end)+vec_ph.f2(1))/2;
% y1 = (equat.f2(end)+equat.f2(1))/2;
% % x = [.2 .6];
% % y = [.1 .5];
% % annotation('textarrow',x,y)
% text(40,2.6,'$\leftarrow\,I_p$','Interpreter', 'Latex','FontSize', 14) 
% 
% %secundario
% equat.f1 = ones(1,length(vec_ph.f1)).*fiLrm.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f1);
% equat.f2 = ones(1,length(vec_ph.f2)).*fiLrm.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,vec_ph.f2);
% plot(rad2deg(vec_ph.f1),n_num*equat.f1, '-.','Color',jetcustom(1,:),'LineWidth',1.5)
% plot(rad2deg(vec_ph.f2),n_num*equat.f2, '-.','Color',jetcustom(1,:),'LineWidth',1.5)
% 
% equat.f1 = ones(1,length(vec_ph.f1)).*fiLrm.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f1);
% equat.f2 = ones(1,length(vec_ph.f2)).*fiLrm.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,vec_ph.f2);
% plot(rad2deg(vec_ph.f1),n_num*equat.f1, '-.','Color',jetcustom(2,:),'LineWidth',1.5)
% plot(rad2deg(vec_ph.f2),n_num*equat.f2, '-.','Color',jetcustom(2,:),'LineWidth',1.5)
% 
% equat.f1 = ones(1,length(vec_ph.f1)).*fiLrm.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f1);
% equat.f2 = ones(1,length(vec_ph.f2)).*fiLrm.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,vec_ph.f2);
% plot(rad2deg(vec_ph.f1),n_num*equat.f1, '-.','Color',jetcustom(3,:),'LineWidth',1.5)
% plot(rad2deg(vec_ph.f2),n_num*equat.f2, '-.','Color',jetcustom(3,:),'LineWidth',1.5)
% 
% % x11 = (vec_ph.f2(end)+vec_ph.f2(1))/2;
% % y11 = (equat.f2(end)+equat.f2(1))/2;
% text(-5,4,'$I_s\,\rightarrow$','Interpreter', 'Latex','FontSize', 14) 
% 
% hold off
% legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10\,$mH'},'Location','best','FontSize', 14,'NumColumns',1);
% 
% xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
% grid on
% grid minor
% xlabel('$\phi\,[^{\circ}]$')
% ylabel('$i_{rms}\,$[A]')
% set(gca, 'FontSize', 20)
% 

























% 
% %calcula traço em que Isw igual a zero
% 
% d_limit_p.f1 = simplify(solve(Ip.f1 == 0, d));
% d_limit_p.f2 = simplify(solve(Ip.f2 == 0, d));
% fd_limit_p.f1 = matlabFunction(d_limit_p.f1, 'vars', {L2,M,phi});
% fd_limit_p.f2 = matlabFunction(d_limit_p.f2, 'vars', {L2,M,phi});
% 
% d_limit_s.f1 = simplify(solve(Is.f1 == 0, d));
% d_limit_s.f2 = simplify(solve(Is.f2 == 0, d));
% fd_limit_s.f1 = matlabFunction(d_limit_s.f1, 'vars', {L1,Ldab,M,phi});
% fd_limit_s.f2 = matlabFunction(d_limit_s.f2, 'vars', {L1,Ldab,M,phi});
% 
% 
% %equacao_limite
% d_limit_p.f1 = Ip.f1 == 0;
% d_limit_p.f2 = Ip.f1 == 0;
% 
% pin.f1 = solve(iME.f1*Vo == Po, Po);
% pin.f2 = solve(iME.f2*Vo == Po, Po);
% 
% fpin.f1 = matlabFunction(pin.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% fpin.f2 = matlabFunction(pin.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% 
% 
% 
% eqp.f1 = solve(subs(pin.f1,phi,solve(d_limit_p.f1 == d,phi)),Po);
% eqp.f2 = solve(subs(pin.f2,phi,solve(d_limit_p.f2 == d,phi)),Po);
% 
% 
% eq_p.f1 = solve(Ip.f1 == 0, d);
% eq_s.f1 = solve(Is.f1 == 0, d);
% eq_pot.f1 = solve(iME.f1*Vo == Po, Po);
% fIp.f1 = matlabFunction(eq_p.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% fIs.f1 = matlabFunction(eq_s.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% fPot.f1 = matlabFunction(eq_pot.f1, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% 
% eq_p.f2 = solve(Ip.f2 == 0, d);
% eq_s.f2 = solve(Is.f2 == 0, d);
% eq_pot.f2 = solve(iME.f2*Vo == Po, Po);
% fIp.f2 = matlabFunction(eq_p.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% fIs.f2 = matlabFunction(eq_s.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% fPot.f2 = matlabFunction(eq_pot.f2, 'vars', {L1,L2,Ldab,M,Vi,d,fs,phi});
% 
% 
% 
% %limite pelo primario
% lim_p.f1 = ones(1,length(ph_n.f1)).*fIp.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,ph_n.f1);
% lim_s.f1 = ones(1,length(ph_n.f1)).*fIs.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,ph_n.f1);
% lim_pot.f1 = ones(1,length(ph_n.f1)).*fPot.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,ph_n.f1);
% lim_p.f2 = ones(1,length(ph_n.f2)).*fIp.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,ph_n.f2);
% lim_s.f2 = ones(1,length(ph_n.f2)).*fIs.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,ph_n.f2);
% lim_pot.f2 = ones(1,length(ph_n.f2)).*fPot.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,ph_n.f2);
% 
% 
% 
% %%
% 
% 
% 
% figure
% cmap = f_create_cmap(3, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% hold on
% fplot(phi*180/pi,fd_limit_p.f1(L2_num(1),M_num(1),phi),intervalo.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
% fplot(phi*180/pi,fd_limit_p.f1(L2_num(2),M_num(2),phi),intervalo.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
% fplot(phi*180/pi,fd_limit_p.f1(L2_num(3),M_num(3),phi),intervalo.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
% 
% fplot(phi*180/pi,fd_limit_p.f2(L2_num(1),M_num(1),phi),intervalo.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% fplot(phi*180/pi,fd_limit_p.f2(L2_num(2),M_num(2),phi),intervalo.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
% fplot(phi*180/pi,fd_limit_p.f2(L2_num(3),M_num(3),phi),intervalo.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
% 
% x = intervalo.f1*180/pi;
% y = [fd_limit_s.f1(L1_num(1),Ldab_num,M_num(1),phi), fd_limit_s.f1(L1_num(1),Ldab_num,M_num(1),phi)];
% line(x, y,'Color',jetcustom(1,:),'LineWidth',1.5);
% fplot(phi*180/pi,fd_limit_s.f2(L1_num(1),Ldab_num,M_num(1),phi),intervalo.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% x = intervalo.f1*180/pi;
% y = [fd_limit_s.f1(L1_num(2),Ldab_num,M_num(2),phi), fd_limit_s.f1(L1_num(2),Ldab_num,M_num(2),phi)];
% line(x, y,'Color',jetcustom(2,:),'LineWidth',1.5);
% fplot(phi*180/pi,fd_limit_s.f2(L1_num(2),Ldab_num,M_num(2),phi),intervalo.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
% x = intervalo.f1*180/pi;
% y = [fd_limit_s.f1(L1_num(3),Ldab_num,M_num(3),phi), fd_limit_s.f1(L1_num(3),Ldab_num,M_num(3),phi)];
% line(x, y,'Color',jetcustom(3,:),'LineWidth',1.5);
% fplot(phi*180/pi,fd_limit_s.f2(L1_num(3),Ldab_num,M_num(3),phi),intervalo.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
% hold off
% grid on
% grid minor
% legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
% xlabel('$\phi\,[^{\circ}]$')
% ylabel('$d$')
% set(gca, 'FontSize', 20)
% ylim([0 3])
% xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
% 
% %% ZVS d por Po
% 
% 
% %%
% figure
% hold on
% fplot(phi*180/pi,fpin.f1(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi),intervalo.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
% fplot(phi*180/pi,fpin.f1(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,phi),intervalo.f1,'Color',jetcustom(2,:),'LineWidth',1.5)
% fplot(phi*180/pi,fpin.f1(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,phi),intervalo.f1,'Color',jetcustom(3,:),'LineWidth',1.5)
% 
% fplot(phi*180/pi,fpin.f2(L1_num(1),L2_num(1),Ldab_num,M_num(1),Vi_num,d_num,fs_num,phi),intervalo.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% fplot(phi*180/pi,fpin.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,phi),intervalo.f2,'Color',jetcustom(2,:),'LineWidth',1.5)
% fplot(phi*180/pi,fpin.f2(L1_num(3),L2_num(3),Ldab_num,M_num(3),Vi_num,d_num,fs_num,phi),intervalo.f2,'Color',jetcustom(3,:),'LineWidth',1.5)
% 
% hold off
% xlim([min(intervalo.f1) max(intervalo.f2)]*180/pi)
% grid on
% grid minor
% legend({'$L_m = 0.7\,$mH','$L_m = 1.4\,$mH','$L_m = 10.0\,$mH'},'Location','best','FontSize', 14)
% set(gca, 'FontSize', 20)
% %%
% 
% %%
% figure
% hold on
% fplot(fPot.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,phi),fIp.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,phi), intervalo.f2)
% fplot(fPot.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,phi),fIs.f2(L1_num(2),L2_num(2),Ldab_num,M_num(2),Vi_num,d_num,fs_num,phi), intervalo.f2)
% hold off
% ylim([0 3])
% %%
% 
% 
% 
% 
% 
% %%
% figure
% hold on
% plot(lim_pot.f1,lim_p.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(lim_pot.f2,lim_p.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(lim_pot.f1,lim_s.f1,'Color',jetcustom(1,:),'LineWidth',1.5)
% plot(lim_pot.f2,lim_s.f2,'Color',jetcustom(1,:),'LineWidth',1.5)
% hold off
% ylim([0 3])





































































%%
% color1 = [0.045 0.245 0.745]; % blue
% color2 = [0.635 0.635 0.635]; % gray
% 
% cmap = f_create_cmap(2, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% M1 = [1; 0; 0; 0; 0; 0];
% M2 = [0; 0; 0; 1; 0; 0];
% is1 = ilL(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)*M1;
% is2 = ilL(Ldab_num, L1_num(1), L2_num(1), M_num(1), phi_num, fs_num, Vi_num, d_num)*M2;
% figure
% hold on
% plot(tf(phi_num,fs_num), is1(:,1)','Color',jetcustom(1,:),'LineWidth',1.5)
% plot(tf(phi_num,fs_num), is2(:,1)','Color',jetcustom(2,:),'LineWidth',1.5)
% hold off
% grid on
% grid minor
% legend({'i_a','i_A'},'Location','best','FontSize', 14)
% xlim([0 1/fs_num])
% hold off
% grid on
% grid minor
% xlim([0 1/fs_num])





























































% 
% clear; close all; clc;
% 
% addpath('utils')
% addpath('utils_transf')
% syms fs Vi d Po dt t real positive
% syms a1 b1 z1 a2 b2 z2 real
% syms xa1 xb1 xz1 xa2 xb2 xz2 real
% syms phi real
% syms Ld1 Ld2 n Lm real positive
% 
% Vo = Vi*d;
% Ts = 1/fs;
% 
% % Parametros para validacao numerica
% Vi_num = 400;
% d_num = 1;
% fs_num = 100e3;
% Ldab_num = 61e-6;
% Ld1_num = 10e-6;
% Ld2_num = 6e-6;
% Lm_num = 700e-6;
% phi_num = deg2rad(50);  %[MUDAR]
% n_num = 5/9;
% Po_num = 2000;
% M_num = Lm_num*n_num;
% L1_num = Ld1_num + Lm_num;
% L2_num = Ld2_num + n_num*n_num*Lm_num;
% k_num = M_num/sqrt(L1_num*L2_num);
% 
% syms L1 L2 Ldab M real positive
% 
% x = sym('x_', [6,1], 'real');
% dx = sym('dx_', [6,1], 'real');
% eq = sym('eq_', [4,1], 'real');
% u = sym('u_', [6,1], 'real');
% s = sym('s_', [6,1], 'real');
% 
% u(1:3) = s(1:3)*Vi;
% u(4:6) = s(4:6)*Vo;
% 
% %corrente no dab
% ild = x(1:2);
% ild = [ild; -sum(ild)];
% dild = dx(1:2);
% dild = [dild; -sum(dild)];
% 
% 
% %corrente no primario do trafo
% il = (ild - [ild(2:3); ild(1)])/3;
% dil = (dild - [dild(2:3); dild(1)])/3;
% 
% %corrente no secundario do trafo
% iL = x(4:5);
% iL = [iL; -sum(iL)];
% diL = dx(4:5);
% diL = [diL; -sum(diL)];
% 
% %definicao das tensoes do indutor acoplado
% vP = + L1*dil + M*diL;
% vS = + M*dil + L2*diL;
% 
% %definicao das tensoes nos indutores do dab
% vLdab = Ldab*dild;
% 
% Tdiff = [1 -1 0; 0 1 -1; -1 0 1];
% 
% m_p = Tdiff*u(1:3) == Tdiff*vLdab + vP; %malha primario
% m_s = Tdiff*u(4:6) == vS; %malha secundario
% 
% eq(1:2) = m_p(1:2);
% eq(3:4) = m_s(1:2);
% 
% %%
% dxs = struct2array(solve(eq, dx,'Real',true,'IgnoreAnalyticConstraints',true)).';
% 
% dxs(3) = -sum(dxs(1) + dxs(2));
% dxs(6) = -sum(dxs(4) + dxs(5));
% dxs = simplify(dxs);
% Ms = equationsToMatrix(dxs, [x;s]);
% 
% As = Ms(:,1:6);
% Bs = Ms(:,7:12);
% 
% %% converter alpha beta
% Tcl = 1/3*[2 -1 -1; 0 sqrt(2) -sqrt(2); 1 1 1]; %Transformada de clark
% Tclm = Tcl(1:2,:); %ignorar nivel zero
% Tclxz = kron(eye(2), Tcl);
% Tclx = kron(eye(2), Tclm);
% 
% Axx = simplify(Tclxz*As*inv(Tclxz));
% Bxx = simplify(Tclxz*Bs*inv(Tclxz)); %a pinv deixa uns numeros feios
% Ax = [Axx(:,1:2),Axx(:,4:5)];
% Bx = [Bxx(:,1:2),Bxx(:,4:5)];
% A = [Ax(1:2,:);Ax(4:5,:)];
% B = [Bx(1:2,:);Bx(4:5,:)];
% 
% 
% %% equivalent circuit
% 
% % syms dxa1 dxb1 dxa2 dxb2 real 
% % syms xa1 xb1 xa2 xb2 real 
% % syms ua1 ub1 ua2 ub2 real 
% % 
% % dxc = [dxa1 dxb1 n*dxa2 n*dxb2]';
% % xc = [xa1 xb1 n*xa2 n*xb2]';
% % uc = [ua1 ub1 ua2/n ub2/n]';
% % 
% % eq_cl = dxc == A*xc+B*uc;
% % eq_cl = simplify(eq_cl);
% % 
% % solve(eq_cl, dxc,'Real',true,'IgnoreAnalyticConstraints',true)
% % eq_c = solve(eq_cl, [ua1 ub1 ua2 ub2],'Real',true,'IgnoreAnalyticConstraints',true)
% 
% 
% 
% %%
% 
% [sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi);
% sf = sf+0.5;
% scl = kron(eye(2), Tclm)*sf; %estados de comutacao, a1b1c1-a2b2c2 to alpha1beta1-alpha2beta2
% 
% ts = simplify(ang)*Ts/(2*pi);
% tf = matlabFunction(ts, 'vars', {phi, fs}); %criar funcao para definir os tempos
% 
% %%
% 
% figure
% subplot(2,1,1)
% stairs(tf(phi_num, fs_num) * 1e6, ([sf sf(:,1)]'+[0 -2 -4 -6 -8 -10]), 'LineWidth', 2)
% xlim([0 1/fs_num]*1e6)
% grid on
% subplot(2,1,2)
% stairs(tf(phi_num, fs_num) * 1e6, ([scl scl(:,1)]'+[0 -2 -4 -6]), 'LineWidth', 2)
% xlim([0 1/fs_num]*1e6)
% grid on
% 
% %% Obter valores de regime permanente
% 
% g = simplify(expm(A*dt));
% h = B*dt;
% xcl = sym('x_', [4,length(ts)], 'real');
% x0 = sym('x0_',[4,1], 'real');
% 
% xcl(:,1) = x0;
% for i=1:length(ts)-1
%     dts = (ts(i+1)-ts(i));
%     xcl(:,i+1) = simplify(subs(g, dt, dts)*xcl(:,i) + subs(h, dt, dts)*scl(:,i));
% end
% 
% x0x = struct2array(solve(xcl(:,1) == -xcl(:,7), x0)).';
% x0s = simplify(pinv(Tclx)*subs(xcl, x0, x0x));
% iSw = matlabFunction(x0s.', 'vars', {Ldab, L1, L2, M, phi, fs, Vi, d});
% 
% derivadas = pinv(Tclx)*B*scl; %derivadas dos estados, indutor e trafo secundario
% 
% %%
% color1 = [0.045 0.245 0.745]; % blue
% color2 = [0.635 0.635 0.635]; % gray
% 
% cmap = f_create_cmap(2, color2, color1);
% colormap(cmap)
% jetcustom = cmap;
% M1 = [1; 0; 0; 0; 0; 0];
% M2 = [0; 0; 0; 1; 0; 0];
% is1 = iSw(Ldab_num, L1_num, L2_num, M_num, phi_num, fs_num, Vi_num, d_num)*M1;
% is2 = iSw(Ldab_num, L1_num, L2_num, M_num, phi_num, fs_num, Vi_num, d_num)*M2;
% figure
% hold on
% plot(tf(phi_num,fs_num), is1(:,1)','Color',jetcustom(1,:),'LineWidth',1.5)
% plot(tf(phi_num,fs_num), is2(:,1)','Color',jetcustom(2,:),'LineWidth',1.5)
% hold off
% grid on
% grid minor
% legend({'i_a','i_A'},'Location','best','FontSize', 14)
% xlim([0 1/fs_num])
% hold off
% grid on
% grid minor
% xlim([0 1/fs_num])