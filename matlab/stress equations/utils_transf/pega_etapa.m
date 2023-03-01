function [etapas] = pega_etapa(sf_p,target)
    etapas = [];
    count = 0;
    for i=1:length(sf_p)
        if (sf_p(:,i) == target)
            count = count + 1;
            etapas(count) = i;
        end
    end
end