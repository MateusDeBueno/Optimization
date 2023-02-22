clear; close all; clc;

load('f_YY_menor60.mat')
load('f_YY_harmonic_menor60.mat')


Vp = 400;
L = 67e-6;
n = 5/9;
d = 1;
fs = 100e3;
phi = 50*pi/180;
k = 1;

f_YY(Vp,L,n,d,fs,phi,k)

% if phi<pi/3
%     M = f_YY_menor60(Vp,L,n,d,fs,phi,k);
%     M_harmonic = f_YY_harmonic_menor60(Vp,L,n,d,fs,phi,k);
% else
%     M = f_YY_maior60(Vp,L,n,d,fs,phi,k);
%     M_harmonic = f_YY_harmonic_maior60(Vp,L,n,d,fs,phi,k);
% end


