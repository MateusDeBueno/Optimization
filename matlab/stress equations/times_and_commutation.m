function [sf_p, sf_s, sf, ang, sec_switch] = times_and_commutation(phi_num,pi)   

    Pi = 3.141592653589793238462643383279;
    if (-Pi <= phi_num && phi_num <= -2*Pi/3)
        situacao = 0;
    elseif (-2*Pi/3 < phi_num && phi_num <= -Pi/3)
        situacao = 1;
    elseif (-Pi/3 < phi_num && phi_num <= 0)
        situacao = 2;
    elseif (0 < phi_num && phi_num <= Pi/3)
        situacao = 3;
    elseif (Pi/3 < phi_num && phi_num <= 2*Pi/3)
        situacao = 4;
    elseif (2*Pi/3 < phi_num && phi_num <= Pi)
        situacao = 5;
    else
        situacao = 6;
    end

    syms phi real
    if (situacao==0) %-180 ate -120
        shift_s = -5;
        shift_ang = 3;
        sec_switch = 8;
    elseif (situacao==1) %-120 ate -60
        shift_s = -3;
        shift_ang = 2;
        sec_switch = 10;
    elseif (situacao==2) %-60 ate 0
        shift_s = -1;
        shift_ang = 1;
        sec_switch = 12;
    elseif (situacao==3) %0 ate 60
        shift_s = 1;
        shift_ang = 0;
        sec_switch = 2;
    elseif (situacao==4) %60 ate 120
        shift_s = 3;
        shift_ang = -1;
        sec_switch = 4;
    elseif (situacao==5) %120 ate 180
        shift_s = 5;
        shift_ang = -2;
        sec_switch = 6;
    end
    
    sf_base = [1 0 1; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1]';
    sf_p = kron(sf_base, [1, 1]); %funcao de comutacao do primario
    sf_s = circshift(sf_p,[0,shift_s]); %funcao de comutacao do secundario
    sf = -0.5+[sf_p; sf_s]; %a1b1c1 to a2b2c2
        
    ang_base = 0:pi/3:2*pi;
    ang_base_s = (ang_base + shift_ang*pi/3)+phi;
    ang = sym('ang_', [2*length(ang_base)-1,1], 'real');
    
    for ii = 1:numel(ang_base)-1
        ang(2*ii-[1 0]) = [ang_base(ii),ang_base_s(ii)];
    end
    ang(end) = 2*pi;
    ang = ang';
end