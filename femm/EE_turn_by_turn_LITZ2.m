clc; clear; close all;

format short eng

addpath('C:\femm42\mfiles');
addpath('utils');
addpath('data');

% dados = [];


awg_data = readtable('awg_table.txt');
awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32

%parametros gerais
L_goal = 65e-6; %indutancia desejada
Imax = 12; %max current, used for minimum turns
Bsat = 0.25; %usado para numero minimo de espiras
f = 100e3; %frequency
A_min = 1; %mm quadrados de cobre, usado para calcular numero d strands
A_max = 4;
Pcore_max = 100;
strands_max_total = 50;

w = 2*pi*f;
Ts = 1/f;
d = 1;
phi = 70*pi/180;
n = 15/27;
t_db = 277e-9;
T__DB_N = t_db/Ts;
Vp = 400;
IN = Vp/(w*L_goal);



%parametros do ferrite
alpha = 1.22;
beta = 2.848;
kc = 78.907;
ki = 4.539113451;
ur = 2100;


sbw = 1.02;
i = 1;


for core_code=1:3
    %core geometry
%     core_code = 3;
    core_geometry = f_get_core(core_code);
    
    core_name = core_geometry.core_name;
    G_carretel = core_geometry.G_carretel;
    G = core_geometry.G;
    E = core_geometry.E;
    C = core_geometry.C;
    B = core_geometry.B;
    A = core_geometry.A;
    D = core_geometry.D;
    Ve = core_geometry.Ve;
    Ae = core_geometry.Ae;
    Ac = core_geometry.Ac;
    Rth = core_geometry.Rth;
    s = 1; %carretel espessura
    core_geometry.s = s;

    for awg_chosen=20:32 %chosen wire
        d_litz = awg_wires(awg_chosen);

        A_litz = pi*(d_litz/2)^2;

        strands_min = ceil(A_min/A_litz);
        strands_max = floor(A_max/A_litz);

        if (strands_max > strands_max_total)
            strands_max = strands_max_total;
        end

        for strands=strands_min:strands_max

            d_target = f_get_d_target(strands,d_litz,sbw);


            N_by_layer_max = floor(G_carretel/d_target);
            layer_max = floor((0.5*(B-E)-s)/d_target);

            Nmin = ceil(Imax*L_goal/(Bsat*Ac));
            Nmax = N_by_layer_max*layer_max;

%             core_geometry.core_name = core_name;
%             core_geometry.G_carretel = G_carretel;core_geometry.G = G;core_geometry.E = E;core_geometry.C = C;core_geometry.B = B;core_geometry.A = A;core_geometry.D = D;core_geometry.s = s;
            wire_geometry.d_litz = d_litz;wire_geometry.strands = strands;wire_geometry.sbw = sbw; wire_geometry.awg_chosen = awg_chosen;
            
            if (Nmin < N_by_layer_max)
                Nmin = N_by_layer_max;
            end
            
%             for N=Nmin:N_by_layer_max:Nmax %somente layers completas
                
                N = 15;
                
                wire_geometry.N = N;
                
                %GET DATA
                output = f_trafo_YY(d,phi,n,T__DB_N,alpha,Vp,Ae,N,f);
                Irsm = output(1,10)*IN;
                I_L_max = output(1,end-1)*IN;
                integrais_L = output(1,end);

                %CORE LOSSES
                deltaB = 2*L_goal*I_L_max/(N*Ae);
                watts_per_m3 = ki/(Ts/2)*abs(deltaB)^(beta-alpha)*integrais_L;
                Pcore = Ve*watts_per_m3;

                %GAT CALCULATION
                g = f_get_gap(core_geometry,N,ur,L_goal);

                if (g < 30.5 && g > 0.025 && Pcore < Pcore_max)

                    vector_i(i) = i;
                    vector_d_litz(i) = d_litz;
                    vector_strands(i) = strands;
                    vector_sbw(i) = sbw;
                    vector_awg_chosen(i) = awg_chosen;
                    vector_N(i) = N;
                    n_wires_to_femm = N*strands;
                    vector_n_wires_to_femm(i) = n_wires_to_femm;
                    vector_core_name(i) = core_name;
                    vector_G_carretel(i) = G_carretel;
                    vector_G(i) = G;
                    vector_E(i) = E;
                    vector_C(i) = C;
                    vector_B(i) = B;
                    vector_A(i) = A;
                    vector_D(i) = D;
                    vector_Ve(i) = Ve;
                    vector_Ae(i) = Ae;
                    vector_Ac(i) = Ac;
                    vector_Rth(i) = Rth;
                    vector_s(i) = s;
                    vector_g(i) = g;
                    vector_Pcore(i) = Pcore;
                    
                    
%                     wire_geometry2.i(i) = i;
%                     
%                     wire_geometry2.d_litz(i) = d_litz;
%                     wire_geometry2.strands(i) = strands;
%                     wire_geometry2.sbw(i) = sbw;
%                     wire_geometry2.awg_chosen(i) = awg_chosen;
%                     wire_geometry2.N(i) = N;
%                     n_wires_to_femm = N*strands;
%                     wire_geometry2.n_wires_to_femm(i) = n_wires_to_femm;
%                     
%                     
%                     core_geometry2.core_name(i) = core_name;
%                     core_geometry2.G_carretel(i) = G_carretel;
%                     core_geometry2.G(i) = G;
%                     core_geometry2.E(i) = E;
%                     core_geometry2.C(i) = C;
%                     core_geometry2.B(i) = B;
%                     core_geometry2.A(i) = A;
%                     core_geometry2.D(i) = D;
%                     core_geometry2.Ve(i) = Ve;
%                     core_geometry2.Ae(i) = Ae;
%                     core_geometry2.Ac(i) = Ac;
%                     core_geometry2.Rth(i) = Rth;
%                     core_geometry2.s(i) = s;
%                     
%                     core_geometry2.g(i) = g;
                    
                    
                    
%                     n_wires_to_femm = strands*N;
%                     wire_geometry.n_wires_to_femm = n_wires_to_femm;
%                     vector_core(i) = core_geometry;
%                     vector_wire(i) = wire_geometry;
%                     vector_g(i) = g;
%                     vector_f(i) = f;
%                     vector_i(i) = i;

                    i = i + 1;
                end
            end
%         end
    end
end

i_total = i

%colocar colunas zeradas tbm
Rtot = 0; 
Pwinding = 0;
Ptot = 0;
Temp = 0;
Ltot = 0;
vector_Rtot = zeros(1,length(vector_Pcore));
vector_Pwinding = zeros(1,length(vector_Pcore));
vector_Ptot = zeros(1,length(vector_Pcore));
vector_Temp = zeros(1,length(vector_Pcore));
vector_Ltot = zeros(1,length(vector_Pcore));


matrix_inicial = [vector_i;vector_d_litz;vector_strands;vector_sbw;vector_awg_chosen;vector_N;...
    vector_n_wires_to_femm;vector_core_name;vector_G_carretel;vector_G;vector_E;...
    vector_C;vector_B;vector_A;vector_D;vector_Ve;vector_Ae;vector_Ac;vector_Rth;...
    vector_s;vector_g;vector_Pcore;vector_Rtot;vector_Pwinding;vector_Ptot;vector_Temp;vector_Ltot]';
DADOS_inicial2323 = array2table(matrix_inicial);
getname = @(x) inputname(1);





names_inicial = [convertCharsToStrings(getname(i)), convertCharsToStrings(getname(d_litz)),...
    convertCharsToStrings(getname(strands)), convertCharsToStrings(getname(sbw)), convertCharsToStrings(getname(awg_chosen)),...
    convertCharsToStrings(getname(N)), convertCharsToStrings(getname(n_wires_to_femm)),...
    convertCharsToStrings(getname(core_name)), convertCharsToStrings(getname(G_carretel)),...
    convertCharsToStrings(getname(G)), convertCharsToStrings(getname(E)), convertCharsToStrings(getname(C)),...
    convertCharsToStrings(getname(B)), convertCharsToStrings(getname(A)),...
    convertCharsToStrings(getname(D)), convertCharsToStrings(getname(Ve)),...
    convertCharsToStrings(getname(Ae)), convertCharsToStrings(getname(Ac)),...
    convertCharsToStrings(getname(Rth)), convertCharsToStrings(getname(s)),...
    convertCharsToStrings(getname(g)), convertCharsToStrings(getname(Pcore)),...
    convertCharsToStrings(getname(Rtot)),convertCharsToStrings(getname(Pwinding)),...
    convertCharsToStrings(getname(Ptot)),convertCharsToStrings(getname(Temp))...
    convertCharsToStrings(getname(Ltot))];
DADOS_inicial2323.Properties.VariableNames = names_inicial;
DADOS_inicial2323 = sortrows(DADOS_inicial2323,'n_wires_to_femm'); %ordenar por mais facil de simular


% xdgfgfdfdg

try 
    load('dados_iniciais.mat');
    initial_condition = nnz(DADOS_inicial.Ptot)+1;
catch
    initial_condition = 1;
end

initial_condition


[M,N] = size(DADOS_inicial);
for i=initial_condition:M
    
    i
    
    core_geometry.core_name = DADOS_inicial.core_name(i);
    core_geometry.G_carretel = DADOS_inicial.G_carretel(i);
    core_geometry.G = DADOS_inicial.G(i);
    core_geometry.E = DADOS_inicial.E(i);
    core_geometry.C = DADOS_inicial.C(i);
    core_geometry.B = DADOS_inicial.B(i);
    core_geometry.A = DADOS_inicial.A(i);
    core_geometry.D = DADOS_inicial.D(i);
    core_geometry.Ve = DADOS_inicial.Ve(i);
    core_geometry.Ae = DADOS_inicial.Ae(i);
    core_geometry.Ac = DADOS_inicial.Ac(i);
    core_geometry.Rth = DADOS_inicial.Rth(i);
    core_geometry.s = DADOS_inicial.s(i);
    core_geometry.Pcore = DADOS_inicial.Pcore(i);
    
    wire_geometry.d_litz = DADOS_inicial.d_litz(i);
    wire_geometry.strands = DADOS_inicial.strands(i);
    wire_geometry.sbw = DADOS_inicial.sbw(i);
    wire_geometry.awg_chosen = DADOS_inicial.awg_chosen(i);
    wire_geometry.N = DADOS_inicial.N(i);
    wire_geometry.n_wires_to_femm = DADOS_inicial.n_wires_to_femm(i);
    
    g = DADOS_inicial.g(i);
    
    %GET RESISTANCE AND INDUCTANCE
    [Rtot,Ltot,vector_n_layers,vector_position] = run_inductor(core_geometry,wire_geometry,g,f);
                
    %WINDING LOSSES
    Pwinding = Rtot*Irsm^2;
    
    %TOTAL LOSSES
    Ptot = Pwinding + core_geometry.Pcore;
    
    %TEMPERATURE
    Temp = Rth*Ptot;
    
    % SAVING DATA
    DADOS_inicial.Rtot(i) = Rtot;
    DADOS_inicial.Ltot(i) = Ltot;
    DADOS_inicial.Pwinding(i) = Pwinding;
    DADOS_inicial.Ptot(i) = Ptot;
    DADOS_inicial.Temp(i) = Temp;

%     dados = [dados; i core_name N d_litz awg_chosen strands n_wires_to_femm n_layers Ltot g Rtot Pwinding Pcore Ptot Temp];
%     save('dados.mat', 'dados')
    
    save('dados_iniciais.mat', 'DADOS_inicial')

    i = i + 1;
end

calculados = nnz(DADOS_inicial.Ptot);

figure
scatter(DADOS_inicial.N(1:calculados),DADOS_inicial.Ptot(1:calculados),[],DADOS_inicial.strands(1:calculados),'filled')
xlabel('N')
ylabel('Ptot')
grid on
hcb=colorbar;










































% 
% 
% 
%  T = struct2table(struct_try)
% 
% core_geometry332.nome(i) = 2;
% 
% sequence = sort(struct_try.n_wires_to_femm,'ascend');
% sort(struct_try,sequence)
% qqq
% 
% 
% i_total = i
% 
% n_wires_to_femm = 0;
% n_layers = 0;
% Ltot = 0;
% Rtot = 0;
% Pwinding = 0;
% Pcore = 0;
% Ptot = 0;
% Temp = 0;
% getname = @(x) inputname(1);
% name = getname(core_name);
% names = [convertCharsToStrings(getname(i)), convertCharsToStrings(getname(core_name)), convertCharsToStrings(getname(N)), convertCharsToStrings(getname(d_litz)), convertCharsToStrings(getname(awg_chosen)), convertCharsToStrings(getname(strands)), convertCharsToStrings(getname(n_wires_to_femm)), convertCharsToStrings(getname(n_layers)), convertCharsToStrings(getname(Ltot)), convertCharsToStrings(getname(g)), convertCharsToStrings(getname(Rtot)), convertCharsToStrings(getname(Pwinding)), convertCharsToStrings(getname(Pcore)), convertCharsToStrings(getname(Ptot)),convertCharsToStrings(getname(Temp))];
% 
% %tenta pegar os dados ja calculados, se n acha, cria um vazio
% try 
%     load('dados.mat');
%     [numRows,numCols] = size(dados);
%     initial_condition = (numRows+1);
% catch
%     dados = [];
%     initial_condition = 1;
% end
% 
% 
% 
% 
% 
% 
% for i=initial_condition:max(vector_i)
%     i
%     core_geometry = vector_core(i);
%     wire_geometry = vector_wire(i);
%     g = vector_g(i);
%     f = vector_f(i);
%     Pcore = vector_Pcore(i);     
%     
%     d_litz = wire_geometry.d_litz;
%     awg_chosen = wire_geometry.awg_chosen
%     strands = wire_geometry.strands
%     N = wire_geometry.N
%     core_name = core_geometry.core_name
%     
%     n_wires_to_femm = strands*N;
%     
%     %GET RESISTANCE AND INDUCTANCE
%     [Rtot,Ltot,n_layers,position] = run_inductor(core_geometry,wire_geometry,g,f);
%                 
%     %WINDING LOSSES
%     Pwinding = Rtot*Irsm^2;
%     
%     %TOTAL LOSSES
%     Ptot = Pwinding + Pcore;
%     
%     %TEMPERATURE
%     Temp = Rth*Ptot;
%     
%     % SAVING DATA
%     dados = [dados; i core_name N d_litz awg_chosen strands n_wires_to_femm n_layers Ltot g Rtot Pwinding Pcore Ptot Temp];
%     save('dados.mat', 'dados')
% end
% 
% 
% 
% DADOS = array2table(dados);
% DADOS.Properties.VariableNames = names;
% 
% figure
% hold on
% scatter(DADOS.strands,DADOS.Ptot,[],DADOS.Temp,'filled')
% xlabel('N')
% ylabel('Ptot')
% hcb=colorbar;
% grid on