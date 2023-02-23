clear; close all; clc;



Vp = 400;
L = 61e-6;
n = 5/9;
d = 1;
fs = 100e3;

% phi = 50*pi/180;
% [P_med, eficiency, Ptotal] = YY_losses(Vp,L,n,d,fs,phi)


ang_min = 40*pi/180;
ang_max = 70*pi/180;
counter = 0;
for phi=ang_min:(ang_max-ang_min)/40:ang_max
    counter = counter + 1;
    [P_med(counter), efficiency(counter), Ptotal] = YY_losses(Vp,L,n,d,fs,phi);
    vector_phi(counter) = phi;
end


figure
scatter(P_med,efficiency,'filled')
grid on

phi_exp = [67 64 60 50 40];
n_exp = [94.56 94.63 94.23 93.78 92.62];

figure
hold on
scatter(vector_phi*180/pi,efficiency,'filled')
scatter(phi_exp,n_exp/100,'filled')
hold off
grid on
