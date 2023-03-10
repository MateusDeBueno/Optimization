function [x_rms,x_mean] = f_rms_mean(derivadas,x0s,ts,estados,estados_divisao)
    syms t real

    [M,N] = size(derivadas);

    x_rms = zeros(M,1);
    for i=estados
        dt = ts(i+1) - ts(i);
        equation_intervalo = (derivadas(:,i).*ones(M,1)*t+x0s(:,i));
        x_rms = x_rms + int(equation_intervalo.^2,[0, dt]);
    end
    x_rms = x_rms/(ts(estados_divisao(end)+1)-ts(estados_divisao(1)));
    x_rms = simplify(sqrt(x_rms));


    x_mean = zeros(M,1);
    for i=estados
        dt = ts(i+1) - ts(i);
        equation_intervalo = (derivadas(:,i).*ones(M,1)*t+x0s(:,i));
        x_mean = x_mean + int(equation_intervalo,[0, dt]);
    end
    x_mean = x_mean/(ts(estados_divisao(end)+1)-ts(estados_divisao(1)));
    x_mean = simplify(x_mean);
end