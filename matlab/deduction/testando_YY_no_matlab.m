clear; close all; clc;

% % tirar fourier das correntes
% https://www.youtube.com/watch?v=ev0juGwUz78
% https://www.geeksforgeeks.org/implementation-of-fourier-series-up-to-n-harmonics-in-matlab/
% https://users.math.msu.edu/users/gnagy/teaching/11-winter/mth235/lslides/L36-235.pdf

syms Vp d n phi fs L t k a_L b_L k_c_L a_t b_t k_c_t N_L Ac_L N_s Ac_trafo real positive
syms I0 real

variaveis =         [Vp,    L,      n,      d,  fs,     phi,        k];
ponto_de_operacao = [400,   67e-6,  5/9,    1,  100e3,  50*pi/180,  2];

n_etapas=6;

T=1/fs;
Vs = Vp*d;

dVp = sym('dVp_', [n_etapas,1],'real'); %tensao aplicada pelo primario
dVs = sym('dVs_', [n_etapas,1],'real'); %tensao palicada pelo secundario
ang = sym('ang_', [n_etapas,1],'real'); %angulo de cada ponto
dt = sym('dT_', [n_etapas,1], 'real'); %tepo de cada intervalo
I = sym('I_', [n_etapas,1], 'real'); %corrente em cada ponto
exp = sym('exp_', [n_etapas,1], 'real'); %expressao de cada etapa

dVp(1) = Vp/3;      dVp(2) = Vp/3;      dVp(3) = 2*Vp/3;
dVp(4) = 2*Vp/3;    dVp(5) = Vp/3;      dVp(6) = Vp/3;

dVs(1) = -Vs/3/n;   dVs(2) = Vs/3/n;    dVs(3) = Vs/3/n;
dVs(4) = 2*Vs/3/n;  dVs(5) = 2*Vs/3/n;  dVs(6) = Vs/3/n;

ang(1) = 0;         ang(2) = phi;       ang(3) = pi/3;
ang(4) = pi/3+phi;  ang(5) = 2*pi/3;    ang(6) = 2*pi/3+phi;

dV = simplify(dVp-dVs);

dt(1:end-1) = (ang(2:end)-ang(1:end-1))*T/(2*pi);
dt(end) = (pi-ang(end))*T/(2*pi);
dt = simplify(dt);

a = dV/L;
for i = 1:n_etapas
    I(i+1) = I(i) + dt(i)*a(i);
end

eq_I1 = I(1) == -I(7);
I1 = simplify(solve(eq_I1, I(1)));
I = simplify(subs(I, I(1), I1)); %resolvendo corrente d cada ponto

exp = a*t+I(1:end-1); %equacao de cada intervalo

%% calcula das correntes de entrada
I_in_med = (int(exp(3),[0, dt(3)]) + int(exp(4),[0, dt(4)]))/(dt(3)+dt(4));
I_in_med = simplify(I_in_med);
I_in_rms = sqrt((int(exp(3)^2,[0, dt(3)]) + int(exp(4)^2,[0, dt(4)]))/(dt(3)+dt(4)));
I_in_rms = simplify(I_in_rms);

%% calcula das correntes de saida
I_out_med = (int(exp(4)/n,[0, dt(4)]) + int(exp(5)/n,[0, dt(5)]))/(dt(4)+dt(5));
I_out_med = simplify(I_out_med);
I_out_rms = sqrt((int((exp(4)/n)^2,[0, dt(4)]) + int((exp(5)/n)^2,[0, dt(5)]))/(dt(4)+dt(5)));
I_out_rms = simplify(I_out_rms);

%% calcula da corrente rms no indutor
I_L_rms = int(exp(1)^2,[0, dt(1)]);
for i=2:n_etapas
    I_L_rms = I_L_rms + int(exp(i)^2,[0, dt(i)]);
end
I_L_rms = simplify(sqrt(I_L_rms*2/T));

time = sym('time_', [n_etapas,1], 'real'); %tempo de cada ponto
time(1) = dt(1);
for i=2:n_etapas
    time(i) = time(i-1)+dt(i);
end

limits = sym('limits_', [n_etapas,2], 'real'); %limites de integracao
limits(1,1) = 0;
limits(1,2) = time(1);
for i=2:n_etapas
    limits(i,1) = time(i-1);
    limits(i,2) = time(i);
end

part = sym('part_', [n_etapas,1], 'real'); %expressao de cada etapa

vetor_t = t*ones(n_etapas,1);
part(1) = exp(1);
% part(2:end) = subs(exp(2:end),vetor_t(1:end-1),(t-time(1:end-1)));
for i=2:n_etapas
    part(i) = subs(exp(i),vetor_t(i-1),(t-time(i-1)));
end



% 
% sym I_L_ck
% I_L_ck = 0;
% for i=1:n_etapas
%     I_L_ck = I_L_ck + int(part(i)*exp(-1i*2*pi*fs*k*t), limits(i,:));
% end
% I_L_ck = I_L_ck/(T/2);
% 
% 
% fddffd


I_L_a_k = 0;
for i=1:n_etapas
    I_L_a_k = I_L_a_k + int(part(i)*cos(k*pi*t/(T/2)), limits(i,:));
end
I_L_a_k = I_L_a_k/(T/2);

I_L_b_k = 0;
for i=1:n_etapas
    I_L_b_k = I_L_b_k + int(part(i)*sin(k*pi*t/(T/2)), limits(i,:));
end
I_L_b_k = I_L_b_k/(T/2);

I_L_c_k = 2*abs(I_L_a_k-1i*I_L_b_k); %cos e sin to exponential
I_L_c_k_rms = I_L_c_k/sqrt(2); %peak to rms

% verificar plot
IL = piecewise((limits(1,1) <= t & t < limits(1,2)), part(1) ...
              ,(limits(2,1) <= t & t < limits(2,2)), part(2) ...
              ,(limits(3,1) <= t & t < limits(3,2)), part(3) ...
              ,(limits(4,1) <= t & t < limits(4,2)), part(4) ...
              ,(limits(5,1) <= t & t < limits(5,2)), part(5) ...
              ,(limits(6,1) <= t & t < limits(6,2)), part(6) ...
              );
IL=simplify(IL);

figure
fplot(subs(IL,variaveis,ponto_de_operacao), [0 0.5/100000])
grid on
%% calcula da corrente rms no trafo
%em serie com indutor, a perda vai ser calculada com as resistencias
%refletidas para o primario

I_trafo_rms = I_L_rms;
I_trafo_c_k = I_L_c_k;
I_trafo_c_k_rms = I_trafo_c_k/sqrt(2);
%% corrente nas chaves do primario
I_sw_p_rms = int(exp(1)^2,[0, dt(1)]);
for i=2:n_etapas
    I_sw_p_rms = I_sw_p_rms + int(exp(i)^2,[0, dt(i)]);
end
I_sw_p_rms = simplify(sqrt(I_sw_p_rms/T));

%% corrente nas chaves do secundario
I_sw_s_rms = int((exp(1)/n)^2,[0, dt(1)]);
for i=2:n_etapas
    I_sw_s_rms = I_sw_s_rms + int((exp(i)/n)^2,[0, dt(i)]);
end
I_sw_s_rms = simplify(sqrt(I_sw_s_rms/T));

%% perdas magneticas no trafo

Wb_trafo = sum(abs(dVs).*dt);
Wb_trafo = simplify(Wb_trafo);
% vpa(subs(Wb_trafo,variaveis,ponto_de_operacao))

% 
% 
% integrais_trafo = 0;
% for i=1:n_etapas
%     integrais_trafo = integrais_trafo + int(abs(dVs(i)/(N_s*Ac_trafo)).^a_t, [0, dt(i)]);
% end
% integrais_trafo = simplify(integrais_trafo)
% 


%% perdas magneticas no indutor

Wb_indutor = max(I)*L*2;
% vpa(subs(Wb_indutor,variaveis,ponto_de_operacao))

% 
% integrais_L = 0;
% for i=1:n_etapas
%     integrais_L = integrais_L + int(abs(dV(i)/(N_L*Ac_L)).^a_t, [0, dt(i)]);
% end
% integrais_L = simplify(integrais_L)
% 


%% criar funcao dos resultados
M = [I_in_med,
    I_in_rms,
    I_out_med,
    I_out_rms,
    I_L_rms,
    I_trafo_rms,
    I_sw_p_rms,
    I_sw_s_rms,
    Wb_trafo,
    Wb_indutor];

M_harmonic = [I_L_c_k,
              I_L_a_k,
              I_L_b_k];

f_YY_menor60 = matlabFunction(M, 'vars', variaveis);
save('f_YY_menor60.mat')

f_YY_harmonic_menor60 = matlabFunction(M_harmonic, 'vars', variaveis);
save('f_YY_harmonic_menor60.mat')

% vpa(subs(M,variaveis,ponto_de_operacao))

variaveis
ponto_de_operacao

vpa(subs(M_harmonic,variaveis,ponto_de_operacao))
% vpa(subs(2*sqrt(I_L_a_k^2+I_L_b_k^2),variaveis,ponto_de_operacao))




