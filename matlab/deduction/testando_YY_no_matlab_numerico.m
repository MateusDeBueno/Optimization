clear; close all; clc;


% variaveis =         [Vp,    L,      n,      d,  fs,     phi];
% ponto_de_operacao = [400,   67e-6,  5/9,    1,  100e3,  50*pi/180];

load('f_YY_menor60.mat')
load('f_YY_harmonic_menor60.mat')


Vp = 400;
L = 67e-6;
n = 5/9;
d = 1;
fs = 100e3;
phi = 50*pi/180;
T = 1/fs;

t = 0:T/1000:T/2;

k_max = 300;

a = zeros(1,k_max);
b = a;
eq = 0;
for k=1:k_max
    M = f_YY_harmonic_menor60(Vp,L,n,d,fs,phi,k);
    a(k) = M(2);
    b(k) = M(3);
    eq = eq+a(k)*cos(k*pi*t/(T/2))+b(k)*sin(k*pi*t/(T/2));
end

for k=1:k_max
    M = f_YY_harmonic_menor60(Vp,L,n,d,fs,phi,k);
    spec(k) = M(1);
end



figure
bar(1:k_max,spec)
grid on


figure
plot(t,eq)
grid on


