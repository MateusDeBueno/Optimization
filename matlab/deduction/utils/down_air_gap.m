function [p] = down_air_gap(p,g)
    p = p - ones(length(p),2).*[0 g];
end

