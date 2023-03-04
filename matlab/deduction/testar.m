clear; close all; clc;

%% Parametros das chaves
dab.sw.fitoff = load("f_fitted_off.mat");
dab.sw.fiton = load("f_fitted_on.mat");

dab.sw.Rds_on = 80e-3;
%% Dados do conversor
Vi = 400;
d = 1;
fs = 100e3;
Ldab = 61e-6;
Ld1 = 2e-6;
Ld2 = 2e-6;
Lm = 700e-6;
phi = deg2rad(50);
n = 5/9;
Vo = d*Vi;

%% Funcoes do YY
dab.YYmaior60 = load('dabYY_functions_maior60.mat');
dab.YYmenor60 = load('dabYY_functions_menor60.mat');

%% potencia de saida
Po = Vo*dab.YYmenor60.f_Iout_med(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo);

%% perdas totais
Pt = dabYY_loss(Ldab,n,Ld1,Ld2,Lm,phi,fs,Vi,Vo,dab)

%%
rendimento = (Po - Pt)/Po;
