function [M,M_harmonic] = f_YY(Vp_num,L_num,n_num,d_num,fs_num,phi_num,k_num)
        
    load('f_YY_menor60.mat')
    load('f_YY_harmonic_menor60.mat')
    load('f_YY_maior60.mat')
    load('f_YY_harmonic_maior60.mat')
    
    
    
    addpath('utils');
    addpath('data');
    
    
    % wire data
    awg_data = readtable('awg_table.txt');
    awg_wires = awg_data.Var3; %all wires in mm %wires between awg1 and awg32
    awg_chosen = 38;
    strand = 180;
    d_l = awg_wires(awg_chosen)*1e-3; % [m]
    
    
    
    Vp_num = 400;
    L_num = 67e-6;
    n_num = 5/9;
    d_num = 1;
    fs_num = 100e3;
    phi_num = 50*pi/180;
    k_num = 1;
    
    
    % parametros trafo e indutor (numero de voltas e espiras em paralelo)
    NL = 14;
    Np = 9;
    Ns = 5;
    Sp = 2;
    Ss = 3;
    
    
    
    
    
    Vp = double(Vp_num);
    L = double(L_num);
    n = double(n_num);
    d = double(d_num);
    fs = double(fs_num);
    phi = double(phi_num);
    
    
    if (phi_num<pi/3 && phi_num>0)
        M = f_YY_menor60(Vp,L,n,d,fs,phi,1);
    else
        M = f_YY_maior60(Vp,L,n,d,fs,phi,1);
    end
    
    I_in_med = M(1);
    I_in_rms = M(2);
    I_out_med = M(3);
    I_out_rms = M(4);
    I_L_rms = M(5);
    I_trf_rms = M(6);
    I_sw_p_rms = M(7);
    I_sw_s_rms = M(8);
    I_sw_p_on = M(9);
    I_sw_s_on = M(10);
    Wb_trf = M(11);
    Wb_indutor = M(12);
    P_med = I_in_med*Vp;
    
    %% Conduction loss
    
    Pcond = f_cond_loss(I_sw_p_rms,I_sw_s_rms);
    
    
    % [ZVS_p, ZVS_s] = f_switch_loss(I_sw_p_on,I_sw_s_on,fs)
    
    
    
    % dsdsd
    %% Copper losses
    Pind_cu = 0;
    Ptrf_cu = 0;
    k_max = 5;
    for i=1:20
    %     k = i
        if (phi_num<pi/3)
            M_harmonic = f_YY_harmonic_menor60(Vp,L,n,d,fs,phi,i)/sqrt(2);
        else
            M_harmonic = f_YY_harmonic_maior60(Vp,L,n,d,fs,phi,i)/sqrt(2);
        end
    
        ik_L = M_harmonic(1);
        ik_trf = M_harmonic(2);
    
        Rk_ind = f_get_resistance_ind(fs_num*i,NL,d_l);
        Rk_trf = f_get_resistance_trf(fs_num*i,Np,Ns,Sp,Ss,d_l);
    
        Pind_cu = Pind_cu + Rk_ind*ik_L^2;
        Ptrf_cu = Ptrf_cu + Rk_trf*ik_trf^2;
    end
    
    Pcu = 3*Pind_cu + Ptrf_cu; % Copper losses
    
    
    %% calculo eficiencia
    Ptotal = Pcond+Pcu;
    eficiency = (P_med - Ptotal)/P_med


end