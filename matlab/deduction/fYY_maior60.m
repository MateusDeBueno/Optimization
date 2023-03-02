clear; close all; clc;

% % tirar fourier das correntes
% https://www.youtube.com/watch?v=ev0juGwUz78
% https://www.geeksforgeeks.org/implementation-of-fourier-series-up-to-n-harmonics-in-matlab/
% https://users.math.msu.edu/users/gnagy/teaching/11-winter/mth235/lslides/L36-235.pdf
% https://rcub.ac.in/econtent/ug/bsc/6sem/BSC%20Sem%20VI%20Physics%20Fourier%20transform.pdf

% syms Vp d n fs L t k a_L b_L k_c_L N_L Ac_L a_trf b_trf k_c_trf N_trf Ac_trf real positive
% syms I0 phi real

syms Vp d n fs L t k real positive
syms I0 phi real

% assumeAlso(phi < pi/3)
% assumeAlso(phi > 0)
% 
% assumeAlso(N_trf > 1)
% assumeAlso(N_L > 1)




% variaveis =         [Vp,  L,     n,   d, fs,    phi,       k,  a_L,   b_L,   k_c_L, N_L, Ac_L, a_trf, b_trf, k_c_trf, N_trf, Ac_trf];
% ponto_de_operacao = [400, 67e-6, 5/9, 1, 100e3, 50*pi/180, 1,  1.394, 2.248, 2.448, 14,  0.01, 1.585, 2.756, .402,    5,     0.01];

variaveis =         [Vp,    L,      n,      d,  fs,     phi,        k];
c
n_etapas=6;

T=1/fs;
Vs = Vp*d;

dVp = sym('dVp_', [n_etapas,1],'real'); %tensao aplicada pelo primario
dVs = sym('dVs_', [n_etapas,1],'real'); %tensao palicada pelo secundario
ang = sym('ang_', [n_etapas,1],'real'); %angulo de cada ponto
dt = sym('dT_', [n_etapas,1], 'real'); %tepo de cada intervalo
I = sym('I_', [n_etapas,1], 'real'); %corrente em cada ponto
exp = sym('exp_', [n_etapas,1], 'real'); %expressao de cada etapa
time = sym('time_', [n_etapas,1], 'real'); %tempo de cada ponto
part = sym('part_', [n_etapas,1], 'real'); %expressao de cada etapa

dVp(1) = Vp/3;      dVp(2) = Vp/3;      dVp(3) = 2*Vp/3;
dVp(4) = 2*Vp/3;    dVp(5) = Vp/3;      dVp(6) = Vp/3;

dVs(1) = -2*Vs/3;   dVs(2) = -Vs/3;     dVs(3) = -Vs/3;
dVs(4) = Vs/3;      dVs(5) = Vs/3;      dVs(6) = 2*Vs/3;
dVs = dVs/n;

ang(1) = 0;         ang(2) = phi-pi/3;  ang(3) = pi/3;
ang(4) = phi;       ang(5) = 2*pi/3;    ang(6) = pi/3+phi;

dV = simplify(dVp-dVs); %tensao applicada no indutor

dt(1:end-1) = (ang(2:end)-ang(1:end-1))*T/(2*pi); %tempo de cada invervalo
dt(end) = (pi-ang(end))*T/(2*pi);
dt = simplify(dt);

a = dV/L;
for i = 1:n_etapas
    I(i+1) = I(i) + dt(i)*a(i);
end
eq_I1 = I(1) == -I(end);
I1 = simplify(solve(eq_I1, I(1))); %corrente inicial no indutor
I = simplify(subs(I, I(1), I1)); %resolvendo corrente d cada ponto

exp = a*t+I(1:end-1); %equacao de cada intervalo no indutor

%% calcula das correntes de entrada
I_in_med = (int(exp(3),[0, dt(3)]) + int(exp(4),[0, dt(4)]))/(dt(3)+dt(4));
I_in_med = simplify(I_in_med);
I_in_rms = sqrt((int(exp(3)^2,[0, dt(3)]) + int(exp(4)^2,[0, dt(4)]))/(dt(3)+dt(4)));
I_in_rms = simplify(I_in_rms);

%% calcula das correntes de saida
I_out_med = (int(exp(6)/n,[0, dt(6)]) + int(-exp(1)/n,[0, dt(1)]))/(dt(6)+dt(1));
I_out_med = simplify(I_out_med);
I_out_rms = sqrt((int((exp(6)/n)^2,[0, dt(6)]) + int((-exp(1)/n)^2,[0, dt(1)]))/(dt(6)+dt(1)));
I_out_rms = simplify(I_out_rms);

%% calcula da corrente rms no indutor

I_L_rms = int(exp(1)^2,[0, dt(1)]);
for i=2:n_etapas
    I_L_rms = I_L_rms + int(exp(i)^2,[0, dt(i)]);
end
I_L_rms = simplify(sqrt(I_L_rms*2/T));

%% calcula os coeficientes da corrente no indutor

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


vetor_t = t*ones(n_etapas,1);
part(1) = exp(1);
for i=2:n_etapas
    part(i) = subs(exp(i),vetor_t(i-1),(t-time(i-1)));
end

% calculo do segundo ciclo
limits2 = limits+simplify(sum(dt));
part2 = -part;
part2 = subs(part2,t,(t-simplify(sum(dt))));
% 
% I_L_a_k = 0;
% for i=1:n_etapas
%     I_L_a_k = I_L_a_k + int(part(i)*cos(2*pi*k*t/T), limits(i,:));
% end
% for i=1:n_etapas
%     I_L_a_k = I_L_a_k + int(part2(i)*cos(2*pi*k*t/T), limits2(i,:));
% end
% I_L_a_k = simplify(I_L_a_k*2/T);
% 
% I_L_b_k = 0;
% for i=1:n_etapas
%     I_L_b_k = I_L_b_k + int(part(i)*sin(2*pi*k*t/T), limits(i,:));
% end
% for i=1:n_etapas
%     I_L_b_k = I_L_b_k + int(part2(i)*sin(2*pi*k*t/T), limits2(i,:));
% end
% I_L_b_k = simplify(I_L_b_k*2/T);
% 
% I_L_c_k = abs(I_L_a_k-1i*I_L_b_k); %cos e sin to exponential
% I_L_c_k = simplify(I_L_c_k);
% I_L_c_k_rms = I_L_c_k/sqrt(2); %peak to rms


% IL = piecewise((limits(1,1) <= t & t < limits(1,2)), part(1) ...
%               ,(limits(2,1) <= t & t < limits(2,2)), part(2) ...
%               ,(limits(3,1) <= t & t < limits(3,2)), part(3) ...
%               ,(limits(4,1) <= t & t < limits(4,2)), part(4) ...
%               ,(limits(5,1) <= t & t < limits(5,2)), part(5) ...
%               ,(limits(6,1) <= t & t < limits(6,2)), part(6) ...
%               );
% IL=simplify(IL);
% IL2 = piecewise((limits2(1,1) <= t & t < limits2(1,2)), part2(1) ...
%               ,(limits2(2,1) <= t & t < limits2(2,2)), part2(2) ...
%               ,(limits2(3,1) <= t & t < limits2(3,2)), part2(3) ...
%               ,(limits2(4,1) <= t & t < limits2(4,2)), part2(4) ...
%               ,(limits2(5,1) <= t & t < limits2(5,2)), part2(5) ...
%               ,(limits2(6,1) <= t & t < limits2(6,2)), part2(6) ...
%               );
% IL2=simplify(IL2);
% figure
% hold on
% fplot(subs(IL,variaveis,ponto_de_operacao), [0 1/100000])
% fplot(subs(IL2,variaveis,ponto_de_operacao), [0 1/100000])
% hold off
% grid on


%% calcula da corrente rms no trf
I_trf_rms = I_L_rms;

%% calcula os coeficientes da corrente no trf
I_trf_c_k = I_L_c_k;
I_trf_c_k_rms = I_trf_c_k/sqrt(2);

%% corrente nas chaves do primario
I_sw_p_rms = int(exp(1)^2,[0, dt(1)]);
for i=2:n_etapas
    I_sw_p_rms = I_sw_p_rms + int(exp(i)^2,[0, dt(i)]);
end
I_sw_p_rms = simplify(sqrt(I_sw_p_rms/T));

I_sw_p_on = I(1);
%% corrente nas chaves do secundario
I_sw_s_rms = int((exp(1)/n)^2,[0, dt(1)]);
for i=2:n_etapas
    I_sw_s_rms = I_sw_s_rms + int((exp(i)/n)^2,[0, dt(i)]);
end
I_sw_s_rms = simplify(sqrt(I_sw_s_rms/T));

I_sw_s_on = -I(4)/n;
%% perdas magneticas no trf
Wb_trf = sum(abs(dVs).*dt);
Wb_trf = simplify(Wb_trf);

% integrais_trf = 0;
% for i=1:n_etapas
%     integrais_trf = integrais_trf + int(abs(dVs(i)/(N_trf*Ac_trf))^a_trf, [0, dt(i)]);
% end
% integrais_trf = simplify(integrais_trf);

%% perdas magneticas no indutor
Wb_indutor = max(I)*L*2;

% integrais_L = 0;
% for i=1:n_etapas
%     integrais_L = integrais_L + int(abs(dV(i)/(N_L*Ac_L)).^a_L, [0, dt(i)]);
% end
% integrais_L = simplify(integrais_L);

%% criar matrizes dos resultados
M = [I_in_med,
    I_in_rms,
    I_out_med,
    I_out_rms,
    I_L_rms,
    I_trf_rms,
    I_sw_p_rms,
    I_sw_s_rms,
    I_sw_p_on,
    I_sw_s_on,
    Wb_trf,
    Wb_indutor];

M_harmonic = [I_L_c_k,
              I_trf_c_k];

vpa(subs(M,variaveis,ponto_de_operacao))
%% exportar funcao

f_YY_maior60 = matlabFunction(M, 'vars', variaveis);
save('f_YY_maior60.mat', "f_YY_maior60", '-mat')

f_YY_harmonic_maior60 = matlabFunction(M_harmonic, 'vars', variaveis);
save('f_YY_harmonic_maior60.mat', "f_YY_harmonic_maior60", '-mat')


