function [p] = horizontal_move(p,x)
    p = p + (ones(length(p),2).*[x 0]);
end