function output = f_trafo_YY(d,phi,n,T__DB_N,alpha,V__p,A__eL,N__L,f)
 
    %Funcao das equacoes do trafo YY
    minimo = 0;
    meio = pi / 0.3e1;
    maximo = pi / 0.2e1;
    
    %frequencia de operacao (so usa nas integrais)
    omega = f*2*pi;
    
    %zera todos os parametros iniciais
    valido = 0;
    Io_med_norm = 0;
    Io_ef_norm = 0;
    FP_in = 0;
    Iin_ef_norm = 0;
    FP_out = 0;
    Ip_entrada_norm = 0;
    Is_entrada_norm = 0;
    FP_trafo = 0;
    IL_ef_norm = 0;
    Ichave_p_ef_norm = 0;
    Ichave_s_ef_norm = 0;
    VL_ef_norm = 0;
    ZVS_check = 0;
    ZVS_db_check = 0;
    I_p_diode_ef_norm = 0;
    I_s_diode_ef_norm = 0;
    I_p_diode_med_norm = 0;
    I_s_diode_med_norm = 0;
    ZVS_primario = 0;
    ZVS_secundario = 0;
    I_L_max_norm = 0;
    integrais_L = 0;
    
    if (phi >= minimo && phi < meio) %angulo menor
        valido = 1;
        Io_med_norm = 0.1e1 / n * phi * (-(3 * phi) + 0.4e1 * pi) / pi / 0.6e1;
        Io_ef_norm = sqrt(0.3e1) / n ^ 2 * pi ^ (-0.1e1 / 0.2e1) * sqrt((-n + d) ^ 2 * pi ^ 3 + 0.27e2 * phi ^ 2 * n * (d + 0.3e1 * n) * pi - 0.54e2 * n * (d + 0.3e1 / 0.2e1 * n) * phi ^ 3) / 0.27e2;
        FP_in = 0.6e1 * ((-n + d) ^ 2 * pi ^ 3 + 0.81e2 * (d + n / 0.3e1) * phi ^ 2 * d * pi - 0.81e2 * phi ^ 3 * d * (d + 0.2e1 / 0.3e1 * n)) ^ (-0.1e1 / 0.2e1) * sqrt(0.3e1) * (-0.3e1 / 0.4e1 * phi + pi) * d * pi ^ (-0.1e1 / 0.2e1) * phi;
        Iin_ef_norm = sqrt(0.3e1) / n * pi ^ (-0.1e1 / 0.2e1) * sqrt((-n + d) ^ 2 * pi ^ 3 + 0.81e2 * (d + n / 0.3e1) * phi ^ 2 * d * pi - 0.81e2 * phi ^ 3 * d * (d + 0.2e1 / 0.3e1 * n)) / 0.27e2;
        FP_out = 0.6e1 * ((-n + d) ^ 2 * pi ^ 3 + 0.27e2 * phi ^ 2 * n * (d + 0.3e1 * n) * pi - 0.54e2 * n * (d + 0.3e1 / 0.2e1 * n) * phi ^ 3) ^ (-0.1e1 / 0.2e1) * sqrt(0.3e1) * pi ^ (-0.1e1 / 0.2e1) * phi * n * (-0.3e1 / 0.4e1 * phi + pi);
        Ip_entrada_norm = (0.2e1 * pi * d - 0.2e1 * pi * n - 0.3e1 * d * phi) / n / 0.9e1;
        Is_entrada_norm = -(0.2e1 * pi * d - 0.2e1 * pi * n + 0.3e1 * n * phi) / n ^ 2 / 0.9e1;
        FP_trafo = 0.3e1 * (0.5e1 * (-n + d) ^ 2 * pi ^ 3 + 0.54e2 * pi * d * n * phi ^ 2 - 0.27e2 * phi ^ 3 * d * n) ^ (-0.1e1 / 0.2e1) * sqrt(0.6e1) * pi ^ (-0.1e1 / 0.2e1) * (-0.3e1 / 0.4e1 * phi + pi) * phi * n;
        IL_ef_norm = sqrt(0.3e1) / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.5e1 * (-n + d) ^ 2 * pi ^ 3 + 0.54e2 * pi * d * n * phi ^ 2 - 0.27e2 * phi ^ 3 * d * n) / 0.27e2;
        Ichave_p_ef_norm = 0.1e1 / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.5e1 * (-n + d) ^ 2 * pi ^ 3 + 0.54e2 * pi * d * n * phi ^ 2 - 0.27e2 * phi ^ 3 * d * n) * sqrt(0.6e1) / 0.54e2;
        Ichave_s_ef_norm = 0.1e1 / n ^ 2 * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.5e1 * (-n + d) ^ 2 * pi ^ 3 + 0.54e2 * pi * d * n * phi ^ 2 - 0.27e2 * phi ^ 3 * d * n) * sqrt(0.6e1) / 0.54e2;
        VL_ef_norm = sqrt(0.2e1) / n * pi ^ (-0.1e1 / 0.2e1) * sqrt((-n + d) ^ 2 * pi + 0.3e1 * d * n * phi) / 0.3e1;
        I_L_max_norm = max([((2 * d - 2 * n) * pi - (3 * d * phi)) / n / 0.9e1 (0.2e1 * pi * d + (-0.2e1 * pi + (3 * phi)) * n) / n / 0.9e1 ((-n + d) * pi + (3 * d * phi)) / n / 0.9e1 ((-pi + (6 * phi)) * n + pi * d) / n / 0.9e1 ((-d + n) * pi + (6 * d * phi)) / n / 0.9e1 ((pi + (3 * phi)) * n - pi * d) / n / 0.9e1]);
        integrais_L = ((2 ^ alpha + 2) * (pi - (3 * phi)) * (abs(-n + d) ^ alpha) / 0.3e1 + ((abs(-n + 2 * d) ^ alpha + (n + d) ^ alpha + abs(-2 * n + d) ^ alpha) * phi)) * V__p ^ alpha * N__L ^ (-alpha) * A__eL ^ (-alpha) * (n ^ (-alpha)) * (3 ^ (-alpha)) / omega;
        
        if (Ip_entrada_norm<0)
            if ((((-2 * d + 2 * n) * pi + (3 * d * phi)) / pi / (n + d) / 0.6e1 < 0.0e0 || T__DB_N < ((-2 * d + 2 * n) * pi + (3 * d * phi)) / pi / (n + d) / 0.6e1 && 0.0e0 < ((-2 * d + 2 * n) * pi + (3 * d * phi)) / pi / (n + d) / 0.6e1) && 0.0e0 <= (0.2e1 * pi * d - 0.2e1 * pi * n + (3 * n * phi)) / n / 0.9e1 || (((2 * d - 2 * n) * pi + (3 * n * phi)) / (-n + d) / pi / 0.6e1 < 0.0e0 || T__DB_N < ((2 * d - 2 * n) * pi + (3 * d * phi)) / (-n + d) / pi / 0.6e1 && 0.0e0 < ((2 * d - 2 * n) * pi + (3 * n * phi)) / (-n + d) / pi / 0.6e1) && (0.2e1 * pi * d - 0.2e1 * pi * n + (3 * n * phi)) / n / 0.9e1 < 0.0e0)
                ZVS_primario = 1;
            end
        end
        
        if (Is_entrada_norm<0)
            if ((((2 * d - 2 * n) * pi + (3 * n * phi)) / (-n + d) / pi / 0.6e1 < 0.0e0 || T__DB_N < ((2 * d - 2 * n) * pi + (3 * n * phi)) / (-n + d) / pi / 0.6e1 && 0.0e0 < ((2 * d - 2 * n) * pi + (3 * n * phi)) / (-n + d) / pi / 0.6e1) && 0.0e0 <= -(pi * d - pi * n + (3 * d * phi)) / (n ^ 2) / 0.9e1 || (((-n + d) * pi + (3 * d * phi)) / pi / (-2 * n + d) / 0.6e1 < 0.0e0 || T__DB_N < ((2 * d - 3 * n) * pi + (6 * n * phi)) / pi / (-2 * n + d) / 0.6e1 && 0.0e0 < ((-n + d) * pi + (3 * d * phi)) / pi / (-2 * n + d) / 0.6e1) && -(pi * d - pi * n + (3 * d * phi)) / (n ^ 2) / 0.9e1 < 0.0e0)
                ZVS_secundario = 1;
            end
        end
        
        if (Ip_entrada_norm<0 && Is_entrada_norm<0)
            ZVS_check = 1;
            if (ZVS_secundario == 1 && ZVS_primario == 1)
                ZVS_db_check = 1;
            end
        end
        
        I_p_diode_ef_norm = Piecewise(T__DB_N < phi / pi / 0.2e1, 0.1e1 / n * sqrt(0.12e2) * sqrt(T__DB_N * (pi ^ 2 * (n + d) ^ 2 * T__DB_N ^ 2 + pi * (n + d) * ((pi - 0.3e1 / 0.2e1 * phi) * d - pi * n) * T__DB_N + ((pi - 0.3e1 / 0.2e1 * phi) * d - pi * n) ^ 2 / 0.3e1)) / 0.9e1, phi / pi / 0.2e1 <= T__DB_N, 0.1e1 / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.12e2 * (-n + d) ^ 2 * T__DB_N * (T__DB_N ^ 2 - T__DB_N + 0.1e1 / 0.3e1) * pi ^ 3 - 0.18e2 * phi * (-n + d) * T__DB_N * d * (T__DB_N - 0.2e1 / 0.3e1) * pi ^ 2 + 0.9e1 * phi ^ 2 * (T__DB_N * d - 0.2e1 / 0.3e1 * d + 0.2e1 / 0.3e1 * n) * d * pi - 0.3e1 * d * phi ^ 3 * n) / 0.9e1);
        I_s_diode_ef_norm = Piecewise(T__DB_N < (pi / 0.3e1 - phi) / pi / 0.2e1, 0.1e1 / n ^ 2 * sqrt(0.12e2) * sqrt(T__DB_N * ((T__DB_N ^ 2 - T__DB_N + 0.1e1 / 0.3e1) * (-n + d) ^ 2 * pi ^ 2 - 0.3e1 / 0.2e1 * n * phi * (-n + d) * (T__DB_N - 0.2e1 / 0.3e1) * pi + 0.3e1 / 0.4e1 * n ^ 2 * phi ^ 2)) / 0.9e1, (pi / 0.3e1 - phi) / pi / 0.2e1 <= T__DB_N, sqrt(0.2e1) / n ^ 2 * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.216e3 * pi ^ 3 * (-0.2e1 * n + d) ^ 2 * T__DB_N ^ 3 - 0.216e3 * pi ^ 2 * (-0.2e1 * n + d) * ((-0.3e1 / 0.2e1 * pi + 0.3e1 * phi) * n + pi * d) * T__DB_N ^ 2 + 0.72e2 * pi * ((-0.3e1 / 0.2e1 * pi + 0.3e1 * phi) * n + pi * d) ^ 2 * T__DB_N + 0.5e1 * n * ((-0.6e1 / 0.5e1 * pi + 0.9e1 / 0.5e1 * phi) * n + d * (pi + 0.3e1 / 0.5e1 * phi)) * (pi - 0.3e1 * phi) ^ 2) / 0.54e2);
        I_p_diode_med_norm = Piecewise(T__DB_N < phi / pi / 0.2e1, (0.36e2 * (((T__DB_N + 0.1e1 / 0.3e1) * d + (T__DB_N - 0.1e1 / 0.3e1) * n) * pi - d * phi / 0.2e1) ^ 2 * Piecewise(((0.6e1 * T__DB_N + 0.2e1) * d + 0.6e1 * T__DB_N * n - 0.2e1 * n) * pi / 0.3e1 - d * phi <= 0.0e0, -1, 0.0e0 < ((0.6e1 * T__DB_N + 0.2e1) * d + 0.6e1 * T__DB_N * n - 0.2e1 * n) * pi / 0.3e1 - d * phi, 1) - 0.8e1 * (Piecewise((0.2e1 * d - 0.2e1 * n) * pi / 0.3e1 - d * phi < 0.0e0, 0, 1) - 0.1e1 / 0.2e1) * ((-n + d) * pi - 0.3e1 / 0.2e1 * d * phi) ^ 2) / pi / n / (n + d) / 0.108e3, phi / pi / 0.2e1 <= T__DB_N, (0.36e2 * (n + d) * ((T__DB_N - 0.1e1 / 0.3e1) * (-n + d) * pi - d * phi / 0.2e1) ^ 2 * Piecewise(0.2e1 * (T__DB_N - 0.1e1 / 0.3e1) * (-n + d) * pi - d * phi < 0.0e0, 0, 1) - 0.4e1 * (-n + d) * ((-n + d) * pi - 0.3e1 / 0.2e1 * d * phi) ^ 2 * Piecewise((0.2e1 * d - 0.2e1 * n) * pi / 0.3e1 - d * phi < 0.0e0, 0, 1) - 0.4e1 * (n + d) * ((-n + d) * pi + 0.3e1 / 0.2e1 * n * phi) ^ 2 * Piecewise((-0.2e1 * d + 0.2e1 * n) * pi / 0.3e1 - n * phi < 0.0e0, 0, 1) - 0.18e2 * (-n + d) * (-0.2e1 / 0.9e1 * ((-n + d) * pi + 0.3e1 / 0.2e1 * n * phi) ^ 2 * Piecewise((0.2e1 * d - 0.2e1 * n) * pi / 0.3e1 + n * phi < 0.0e0, 0, 1) + (n + d) * ((T__DB_N - 0.2e1 / 0.3e1) * (-n + d) * T__DB_N * pi ^ 2 - ((T__DB_N - 0.2e1 / 0.3e1) * d + 0.2e1 / 0.3e1 * n) * phi * pi + n * phi ^ 2 / 0.2e1))) / n / (n + d) / (-n + d) / pi / 0.54e2);
        I_s_diode_med_norm = Piecewise(T__DB_N < (pi / 0.3e1 - phi) / pi / 0.2e1, (0.36e2 * ((T__DB_N - 0.1e1 / 0.3e1) * (-n + d) * pi - n * phi / 0.2e1) ^ 2 * Piecewise(0.2e1 * (T__DB_N - 0.1e1 / 0.3e1) * (-n + d) * pi - n * phi <= 0.0e0, -1, 0.0e0 < 0.2e1 * (T__DB_N - 0.1e1 / 0.3e1) * (-n + d) * pi - n * phi, 1) - 0.8e1 * ((-n + d) * pi + 0.3e1 / 0.2e1 * n * phi) ^ 2 * (Piecewise((-0.2e1 * d + 0.2e1 * n) * pi / 0.3e1 - n * phi < 0.0e0, 0, 1) - 0.1e1 / 0.2e1)) / pi / n ^ 2 / (-n + d) / 0.108e3, (pi / 0.3e1 - phi) / pi / 0.2e1 <= T__DB_N, (0.72e2 * (((-0.2e1 * T__DB_N + 0.1e1 / 0.2e1) * pi - phi) * n + pi * d * (T__DB_N - 0.1e1 / 0.3e1)) ^ 2 * (-n + d) * Piecewise(0.2e1 * pi * (-0.2e1 * n + d) * T__DB_N - 0.2e1 / 0.3e1 * pi * d + (pi - 0.2e1 * phi) * n < 0.0e0, 0, 1) - 0.8e1 * (-0.2e1 * n + d) * ((-pi + 0.3e1 / 0.2e1 * phi) * n + pi * d) ^ 2 * Piecewise((-0.2e1 * d + 0.2e1 * n) * pi / 0.3e1 - n * phi < 0.0e0, 0, 1) - 0.2e1 * n * ((pi + 0.3e1 * phi) * d - pi * n) ^ 2 * Piecewise((-d + n) * pi / 0.3e1 - d * phi < 0.0e0, 0, 1) - 0.36e2 * (-n + d) * (-0.2e1 * n + d) * (((-0.2e1 * T__DB_N ^ 2 + T__DB_N - 0.1e1 / 0.36e2) * pi ^ 2 + phi * (-0.2e1 * T__DB_N + 0.1e1 / 0.6e1) * pi - phi ^ 2 / 0.4e1) * n + pi ^ 2 * d * T__DB_N * (T__DB_N - 0.2e1 / 0.3e1))) / n ^ 2 / (-n + d) / (-0.2e1 * n + d) / pi / 0.108e3);
    elseif (phi >= meio && phi <= maximo) %angulo maior
        valido = 1;
        Io_med_norm = -0.1e1 / n * (pi ^ 2 - 0.18e2 * pi * phi + 0.18e2 * phi ^ 2) / pi / 0.18e2;
        Io_ef_norm = sqrt(0.3e1) / n ^ 2 * pi ^ (-0.1e1 / 0.2e1) * sqrt((-0.11e2 * pi ^ 3 + 0.81e2 * pi ^ 2 * phi - 0.81e2 * pi * phi ^ 2) * n ^ 2 + 0.9e1 * n * (pi - 0.2e1 * phi) * (pi ^ 2 - 0.6e1 * pi * phi + 0.6e1 * phi ^ 2) * d + pi ^ 3 * d ^ 2) / 0.27e2;
        FP_in = -d * (pi ^ 2 - 0.18e2 * pi * phi + 0.18e2 * phi ^ 2) * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.3e1) * ((-0.11e2 * pi ^ 3 + 0.81e2 * pi ^ 2 * phi - 0.81e2 * pi * phi ^ 2) * d ^ 2 + 0.9e1 * n * (pi - 0.2e1 * phi) * (pi ^ 2 - 0.6e1 * pi * phi + 0.6e1 * phi ^ 2) * d + pi ^ 3 * n ^ 2) ^ (-0.1e1 / 0.2e1) / 0.2e1;
        Iin_ef_norm = 0.1e1 / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.3e1) * sqrt((-0.11e2 * pi ^ 3 + 0.81e2 * pi ^ 2 * phi - 0.81e2 * pi * phi ^ 2) * d ^ 2 + 0.9e1 * n * (pi - 0.2e1 * phi) * (pi ^ 2 - 0.6e1 * pi * phi + 0.6e1 * phi ^ 2) * d + pi ^ 3 * n ^ 2) / 0.27e2;
        FP_out = -n * (pi ^ 2 - 0.18e2 * pi * phi + 0.18e2 * phi ^ 2) * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.3e1) * ((-0.11e2 * pi ^ 3 + 0.81e2 * pi ^ 2 * phi - 0.81e2 * pi * phi ^ 2) * n ^ 2 + 0.9e1 * n * (pi - 0.2e1 * phi) * (pi ^ 2 - 0.6e1 * pi * phi + 0.6e1 * phi ^ 2) * d + pi ^ 3 * d ^ 2) ^ (-0.1e1 / 0.2e1) / 0.2e1;
        Ip_entrada_norm = (0.3e1 * pi * d - 0.2e1 * pi * n - 0.6e1 * d * phi) / n / 0.9e1;
        Is_entrada_norm = -(0.2e1 * pi * d - 0.3e1 * pi * n + 0.6e1 * n * phi) / n ^ 2 / 0.9e1;
        FP_trafo = -(pi ^ 2 - 0.18e2 * pi * phi + 0.18e2 * phi ^ 2) * n * pi ^ (-0.1e1 / 0.2e1) * (0.5e1 * pi ^ 3 * d ^ 2 - 0.9e1 * n * (pi - 0.2e1 * phi) * (pi ^ 2 + 0.3e1 * pi * phi - 0.3e1 * phi ^ 2) * d + 0.5e1 * pi ^ 3 * n ^ 2) ^ (-0.1e1 / 0.2e1) * sqrt(0.6e1) / 0.4e1;
        IL_ef_norm = sqrt(0.3e1) / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.5e1 * pi ^ 3 * d ^ 2 - 0.9e1 * n * (pi - (2 * phi)) * (pi ^ 2 + 0.3e1 * pi * phi - (3 * phi ^ 2)) * d + 0.5e1 * pi ^ 3 * n ^ 2) / 0.27e2;
        Ichave_p_ef_norm = 0.1e1 / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.5e1 * pi ^ 3 * d ^ 2 - 0.9e1 * n * (pi - (2 * phi)) * (pi ^ 2 + 0.3e1 * pi * phi - (3 * phi ^ 2)) * d + 0.5e1 * pi ^ 3 * n ^ 2) * sqrt(0.6e1) / 0.54e2;
        Ichave_s_ef_norm = 0.1e1 / n ^ 2 * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.5e1 * pi ^ 3 * d ^ 2 - 0.9e1 * n * (pi - (2 * phi)) * (pi ^ 2 + 0.3e1 * pi * phi - (3 * phi ^ 2)) * d + 0.5e1 * pi ^ 3 * n ^ 2) * sqrt(0.6e1) / 0.54e2;
        VL_ef_norm = sqrt(0.2e1) / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(pi * d ^ 2 - 0.3e1 * n * (pi - (2 * phi)) * d + pi * n ^ 2) / 0.3e1;
        I_L_max_norm = max([((3 * d - 2 * n) * pi - (6 * d * phi)) / n / 0.9e1 ((-0.3e1 * pi + (3 * phi)) * n + pi * d) / n / 0.9e1 ((3 * d - n) * pi - (3 * d * phi)) / n / 0.9e1 (0.2e1 * pi * d + (-0.3e1 * pi + (6 * phi)) * n) / n / 0.9e1 (n * pi + (3 * d * phi)) / n / 0.9e1 (pi * d + (3 * n * phi)) / n / 0.9e1]);
        integrais_L = 0.2e1 / 0.3e1 * A__eL ^ (-alpha) * V__p ^ alpha * (3 ^ (-alpha)) * (n ^ (-alpha)) * N__L ^ (-alpha) * ((pi - 0.3e1 / 0.2e1 * phi) * (abs(-2 * n + d) ^ alpha) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * (abs(d - n) ^ alpha) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * ((2 * n + d) ^ alpha) + (-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * ((2 * d + n) ^ alpha) + (pi - 0.3e1 / 0.2e1 * phi) * (abs(-n + 2 * d) ^ alpha + (n + d) ^ alpha)) / omega;

        
        if (Ip_entrada_norm<0)
            if ((((-3 * d + 2 * n) * pi + (6 * d * phi)) / (6 * n + 12 * d) / pi < 0.0e0 || T__DB_N < ((-3 * d + 2 * n) * pi + (6 * d * phi)) / (6 * n + 12 * d) / pi && 0.0e0 < ((-3 * d + 2 * n) * pi + (6 * d * phi)) / (6 * n + 12 * d) / pi) && 0.0e0 <= (pi * d - 0.3e1 * pi * n + (3 * n * phi)) / n / 0.9e1 || (((-d + 3 * n) * pi - (3 * n * phi)) / pi / (n + d) / 0.6e1 < 0.0e0 || T__DB_N < ((-2 * d + 2 * n) * pi + (3 * d * phi)) / pi / (n + d) / 0.6e1 && 0.0e0 < ((-d + 3 * n) * pi - (3 * n * phi)) / pi / (n + d) / 0.6e1) && (pi * d - 0.3e1 * pi * n + (3 * n * phi)) / n / 0.9e1 < 0.0e0)
                ZVS_primario = 1;
            end
        end
        
        if (Is_entrada_norm<0)
            if ((((-3 * n + 2 * d) * pi + (6 * n * phi)) / pi / (-2 * n + d) / 0.6e1 < 0.0e0 || T__DB_N < ((-3 * n + 2 * d) * pi + (6 * n * phi)) / pi / (-2 * n + d) / 0.6e1 && 0.0e0 < ((-3 * n + 2 * d) * pi + (6 * n * phi)) / pi / (-2 * n + d) / 0.6e1) && 0.0e0 <= -(pi * n + (3 * d * phi)) / (n ^ 2) / 0.9e1 || ((pi * n + (3 * d * phi)) / pi / (-n + d) / 0.6e1 < 0.0e0 || T__DB_N < ((-n + 2 * d) * pi + (3 * n * phi)) / (-n + d) / pi / 0.6e1 && 0.0e0 < (pi * n + (3 * d * phi)) / pi / (-n + d) / 0.6e1) && -(pi * n + (3 * d * phi)) / (n ^ 2) / 0.9e1 < 0.0e0)
                ZVS_secundario = 1;
            end
        end
        
        if (Ip_entrada_norm<0 && Is_entrada_norm<0)
            ZVS_check = 1;
            if (ZVS_secundario == 1 && ZVS_primario == 1)
                ZVS_db_check = 1;
            end
        end
        
        I_p_diode_ef_norm = Piecewise(T__DB_N < (phi - pi / 0.3e1) / pi / 0.2e1, 0.1e1 / n * sqrt(0.48e2) * sqrt(T__DB_N * ((d + n / 0.2e1) ^ 2 * pi ^ 2 * T__DB_N ^ 2 + 0.3e1 / 0.4e1 * (d + n / 0.2e1) * ((pi - 0.2e1 * phi) * d - 0.2e1 / 0.3e1 * pi * n) * pi * T__DB_N + 0.3e1 / 0.16e2 * ((pi - 0.2e1 * phi) * d - 0.2e1 / 0.3e1 * pi * n) ^ 2)) / 0.9e1, (phi - pi / 0.3e1) / pi / 0.2e1 <= T__DB_N, sqrt(0.2e1) / n * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.216e3 * pi ^ 3 * (n + d) ^ 2 * T__DB_N ^ 3 + 0.216e3 * (n + d) * ((pi - 0.3e1 / 0.2e1 * phi) * d - pi * n) * pi ^ 2 * T__DB_N ^ 2 + 0.72e2 * ((pi - 0.3e1 / 0.2e1 * phi) * d - pi * n) ^ 2 * pi * T__DB_N - 0.6e1 * (pi - 0.3e1 * phi) ^ 2 * d * ((pi - 0.3e1 / 0.2e1 * phi) * d - 0.7e1 / 0.6e1 * n * (pi - 0.3e1 / 0.7e1 * phi))) / 0.54e2);
        I_s_diode_ef_norm = Piecewise(T__DB_N < (0.2e1 / 0.3e1 * pi - phi) / pi / 0.2e1, 0.1e1 / (n ^ 2) * sqrt(0.12e2) * sqrt(T__DB_N * (pi ^ 2 * ((-2 * n + d) ^ 2) * T__DB_N ^ 2 - (-2 * n + d) * ((-0.3e1 / 0.2e1 * pi + 0.3e1 * phi) * n + pi * d) * pi * T__DB_N + ((-0.3e1 / 0.2e1 * pi + 0.3e1 * phi) * n + pi * d) ^ 2 / 0.3e1)) / 0.9e1, (0.2e1 / 0.3e1 * pi - phi) / pi / 0.2e1 <= T__DB_N, sqrt(0.2e1) / (n ^ 2) * pi ^ (-0.1e1 / 0.2e1) * sqrt(0.216e3 * pi ^ 3 * ((-n + d) ^ 2) * T__DB_N ^ 3 - 0.216e3 * (-n + d) * ((-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * n + pi * d) * pi ^ 2 * T__DB_N ^ 2 + 0.72e2 * ((-pi / 0.2e1 + 0.3e1 / 0.2e1 * phi) * n + pi * d) ^ 2 * pi * T__DB_N - 0.16e2 * ((-0.3e1 / 0.4e1 * pi + 0.9e1 / 0.4e1 * phi) * n + d * (pi + 0.3e1 / 0.4e1 * phi)) * n * (pi - 0.3e1 / 0.2e1 * phi) ^ 2) / 0.54e2);
        I_p_diode_med_norm = Piecewise(T__DB_N < (phi - pi / 0.3e1) / pi / 0.2e1, (0.144e3 * (((T__DB_N + 0.1e1 / 0.4e1) * d + (T__DB_N - 0.1e1 / 0.3e1) * n / 0.2e1) * pi - d * phi / 0.2e1) ^ 2 * Piecewise(((0.12e2 * T__DB_N + 0.3e1) * d + 0.6e1 * T__DB_N * n - 0.2e1 * n) * pi / 0.3e1 - 0.2e1 * d * phi <= 0.0e0, -1, 0.0e0 < ((0.12e2 * T__DB_N + 0.3e1) * d + 0.6e1 * T__DB_N * n - 0.2e1 * n) * pi / 0.3e1 - 0.2e1 * d * phi, 1) - 0.18e2 * (Piecewise((0.3e1 * d - 0.2e1 * n) * pi / 0.3e1 - 0.2e1 * d * phi < 0.0e0, 0, 1) - 0.1e1 / 0.2e1) * ((d - 0.2e1 / 0.3e1 * n) * pi - 0.2e1 * d * phi) ^ 2) / pi / n / (n + 0.2e1 * d) / 0.108e3, (phi - pi / 0.3e1) / pi / 0.2e1 <= T__DB_N, (0.144e3 * (((T__DB_N + 0.1e1 / 0.3e1) * pi - phi / 0.2e1) * d + pi * n * (T__DB_N - 0.1e1 / 0.3e1)) ^ 2 * (d + n / 0.2e1) * Piecewise(((0.6e1 * T__DB_N + 0.2e1) * d + 0.6e1 * T__DB_N * n - 0.2e1 * n) * pi / 0.3e1 - d * phi < 0.0e0, 0, 1) - 0.18e2 * ((pi - 0.2e1 * phi) * d - 0.2e1 / 0.3e1 * pi * n) ^ 2 * (n + d) * Piecewise((0.3e1 * d - 0.2e1 * n) * pi / 0.3e1 - 0.2e1 * d * phi < 0.0e0, 0, 1) - 0.2e1 * (pi * d - 0.3e1 * (pi - phi) * n) ^ 2 * d * Piecewise((d - 0.3e1 * n) * pi / 0.3e1 + n * phi < 0.0e0, 0, 1) - 0.72e2 * (n + d) * (((T__DB_N ^ 2 + 0.2e1 / 0.3e1 * T__DB_N - 0.1e1 / 0.36e2) * pi ^ 2 - phi * (T__DB_N - 0.1e1 / 0.6e1) * pi - phi ^ 2 / 0.4e1) * d + pi ^ 2 * n * T__DB_N * (T__DB_N - 0.2e1 / 0.3e1)) * (d + n / 0.2e1)) / (n + d) / n / (d + n / 0.2e1) / pi / 0.216e3);
        I_s_diode_med_norm = Piecewise(T__DB_N < (0.2e1 / 0.3e1 * pi - phi) / pi / 0.2e1, (0.36e2 * (((-0.2e1 * T__DB_N + 0.1e1 / 0.2e1) * n + d * (T__DB_N - 0.1e1 / 0.3e1)) * pi - n * phi) ^ 2 * Piecewise(0.2e1 * pi * (-0.2e1 * n + d) * T__DB_N - 0.2e1 / 0.3e1 * pi * d + (pi - 0.2e1 * phi) * n <= 0.0e0, -1, 0.0e0 < 0.2e1 * pi * (-0.2e1 * n + d) * T__DB_N - 0.2e1 / 0.3e1 * pi * d + (pi - 0.2e1 * phi) * n, 1) - 0.8e1 * ((-0.3e1 / 0.2e1 * n + d) * pi + 0.3e1 * n * phi) ^ 2 * (Piecewise((-0.2e1 * d + 0.3e1 * n) * pi / 0.3e1 - 0.2e1 * n * phi < 0.0e0, 0, 1) - 0.1e1 / 0.2e1)) / pi / n ^ 2 / (-0.2e1 * n + d) / 0.108e3, (0.2e1 / 0.3e1 * pi - phi) / pi / 0.2e1 <= T__DB_N, (0.72e2 * (-0.2e1 * n + d) * (((-T__DB_N + 0.1e1 / 0.6e1) * pi - phi / 0.2e1) * n + pi * d * (T__DB_N - 0.1e1 / 0.3e1)) ^ 2 * Piecewise(((-0.6e1 * T__DB_N + 0.1e1) * n + 0.6e1 * T__DB_N * d - 0.2e1 * d) * pi / 0.3e1 - n * phi < 0.0e0, 0, 1) - 0.36e2 * (0.2e1 / 0.9e1 * ((-0.3e1 / 0.2e1 * pi + 0.3e1 * phi) * n + pi * d) ^ 2 * Piecewise((-0.2e1 * d + 0.3e1 * n) * pi / 0.3e1 - 0.2e1 * n * phi < 0.0e0, 0, 1) + (((-T__DB_N ^ 2 + T__DB_N / 0.3e1 + 0.1e1 / 0.9e1) * pi ^ 2 - phi * (T__DB_N + 0.1e1 / 0.3e1) * pi + phi ^ 2 / 0.4e1) * n + pi ^ 2 * d * T__DB_N * (T__DB_N - 0.2e1 / 0.3e1)) * (-0.2e1 * n + d)) * (-n + d)) / n ^ 2 / (-0.2e1 * n + d) / (-n + d) / pi / 0.108e3);
    end
    
    output = [valido,Io_med_norm,Io_ef_norm,FP_in,Iin_ef_norm,FP_out,Ip_entrada_norm,Is_entrada_norm,FP_trafo,IL_ef_norm,Ichave_p_ef_norm,Ichave_s_ef_norm,VL_ef_norm,ZVS_check,ZVS_db_check,I_p_diode_ef_norm,I_s_diode_ef_norm,I_p_diode_med_norm,I_s_diode_med_norm,ZVS_primario,ZVS_secundario,I_L_max_norm,integrais_L];
    
end