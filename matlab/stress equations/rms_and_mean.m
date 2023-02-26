function [x_rms,x_mean] = rms_and_mean(derivadas,x0s,ts,estados,estados_divisao)
    syms t real
    
    % calcula tempo
    total_time = 0;
    for i=estados_divisao
        dt = ts(i+1) - ts(i);
        total_time = total_time + dt;
    end    
    
    [M,N] = size(derivadas);

    % calcula rms
    x_rms = zeros(M,1);
    for i=estados
        dt = ts(i+1) - ts(i);
        equation_intervalo = (derivadas(:,i).*ones(M,1)*t+x0s(:,i));
        x_rms = x_rms + int(equation_intervalo.^2,[0, dt]);
    end
    x_rms = x_rms/total_time;
    x_rms = simplify(sqrt(x_rms));
    
    % calcula mean
    x_mean = zeros(M,1);
    for i=estados
        dt = ts(i+1) - ts(i);
        equation_intervalo = (derivadas(:,i).*ones(M,1)*t+x0s(:,i));
        x_mean = x_mean + int(equation_intervalo,[0, dt]);
    end
    x_mean = x_mean/total_time;
    x_mean = simplify(x_mean);
end