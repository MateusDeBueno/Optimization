clear; close all; clc;


load('equationss.mat')


%%
Vi = 400;
d = 1;
fs = 100e3;
Ldab = 61e-6;
Ld1 = 2e-6;
Ld2 = 2e-6;
Lm = 700e-6;
phi = deg2rad(-80);  %[MUDAR]
n = 5/9;
Vo = d*Vi;
Po = 3000;
pi = 3.141592653589793;